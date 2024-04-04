
// This code is based on three.js, which comes with the following license:
//
// The MIT License
//
// Copyright © 2010-2024 three.js authors
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
import * as THREE from 'three';

import { GUI } from 'three/addons/libs/lil-gui.module.min.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

let name = 'ToricFresnelToy';

// variables that are not uniforms
let aspectRatioVideoFeedU = 4.0/3.0;
let aspectRatioVideoFeedE = 4.0/3.0;

// Nokia HR20, according to https://www.camerafv5.com/devices/manufacturers/hmd_global/nokia_xr20_ttg_0/
let fovVideoFeedU = 67.3;	// (user-facing) camera
let fovVideoFeedE = 68.3;	// (environment-facing) camera
// the FOV of the screen depends on the user's distance from the screen, of course;
// (randomly) make it similar to the FOV of the cameras (then the camera image fills the screen, which is slightly heartwarming)
let fovScreen = 68;

let cameraLensDistance = 10.0;
let offset = new THREE.Vector3(0, 0, 0.0001);
let deltaPhi = 0.0;

// camera with wide aperture
let apertureRadius = 0.0;
let focusDistance = 1e8;
let noOfRays = 1;

let raytracingSphereRadius = 10000.0;	// what this code does should be independent of this value, but isn't (?)

// "internal" variables
let raytracingSphere;
let raytracingSphereShaderMaterial;

let scene;
let renderer;
let videoFeedU, videoFeedE;	// feeds from user/environment-facing cameras
let camera;
let controls;
let gui;
	
// the status text area
let status = document.createElement('div');
let statusTime;	// the time the last status was posted

// the info text area
let info = document.createElement('div');

// the stored photo
let storedPhoto;
let storedPhotoDescription;
let storedPhotoInfoString;
let showingStoredPhoto = false;

// my Canon EOS450D
const click = new Audio('./click.m4a');

init();
animate();

function init() {
	// create the info element first so that any problems can be communicated
	createStatus();

	scene = new THREE.Scene();
	// scene.background = new THREE.Color( 'skyblue' );
	let windowAspectRatio = window.innerWidth / window.innerHeight;
	camera = new THREE.PerspectiveCamera( fovScreen, windowAspectRatio, 0.5*raytracingSphereRadius, 2*raytracingSphereRadius );
	camera.position.z = cameraLensDistance;
	screenChanged();
	
	renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
	renderer.setPixelRatio(window.devicePixelRatio);
	renderer.setSize( window.innerWidth, window.innerHeight );
	document.body.appendChild( renderer.domElement );
	// document.getElementById('livePhoto').appendChild( renderer.domElement );

	createVideoFeeds();

	addRaytracingSphere();

	// user interface

	addEventListenersEtc();

	addOrbitControls();

	// the controls menu
	createGUI();

	createInfo();
	refreshInfo();
}

function animate() {
	requestAnimationFrame( animate );

	// stats.begin();

	if(!showingStoredPhoto) {
		// update uniforms
		updateUniforms();

		// ensure the raytracing  sphere is centred on  the camera
		raytracingSphere.position.copy(camera.position);

		renderer.render( scene,  camera );
	}

	// stats.end();
}

function updateUniforms() {
	// the tangents for the environment-facing camera video feed
	let tanHalfFovHE, tanHalfFovVE;
	if(aspectRatioVideoFeedE > 1.0) {
		// horizontal orientation
		tanHalfFovHE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0);
		tanHalfFovVE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0)/aspectRatioVideoFeedE;
	} else {
		// vertical orientation
		tanHalfFovHE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0)*aspectRatioVideoFeedE;
		tanHalfFovVE = Math.tan(0.5*fovVideoFeedE*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.halfWidthE.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovHE;
	raytracingSphereShaderMaterial.uniforms.halfHeightE.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovVE;

	// the tangents for the user-facing camera video feed
	let tanHalfFovHU, tanHalfFovVU;
	if(aspectRatioVideoFeedU > 1.0) {
		// horizontal orientation
		tanHalfFovHU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0);
		tanHalfFovVU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0)/aspectRatioVideoFeedU;
	} else {
		// vertical orientation
		tanHalfFovHU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0)*aspectRatioVideoFeedU;
		tanHalfFovVU = Math.tan(0.5*fovVideoFeedU*Math.PI/180.0);
	}
	raytracingSphereShaderMaterial.uniforms.halfWidthU.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovHU;
	raytracingSphereShaderMaterial.uniforms.halfHeightU.value = raytracingSphereShaderMaterial.uniforms.videoDistance.value*tanHalfFovVU;

	// arrange them symmetrically around z=0
	raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.x =  0.5*offset.x;
	raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.y =  0.5*offset.y;
	raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.z =  0.5*offset.z;
	raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.x = -0.5*offset.x;
	raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.y = -0.5*offset.y;
	raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.z = -0.5*offset.z;

	// rotation
	let alpha1 =  0.5*deltaPhi;
	let alpha2 = -0.5*deltaPhi;
	raytracingSphereShaderMaterial.uniforms.cosAlpha1.value = Math.cos(alpha1);
	raytracingSphereShaderMaterial.uniforms.sinAlpha1.value = Math.sin(alpha1);
	raytracingSphereShaderMaterial.uniforms.cosAlpha2.value = Math.cos(alpha2);
	raytracingSphereShaderMaterial.uniforms.sinAlpha2.value = Math.sin(alpha2);

	// create the points on the aperture

	// create basis vectors for the camera's clear aperture
	let viewDirection = new THREE.Vector3();
	let apertureBasisVector1 = new THREE.Vector3();
	let apertureBasisVector2 = new THREE.Vector3();
	camera.getWorldDirection(viewDirection);
	// if(counter < 10) console.log(`viewDirection = (${viewDirection.x.toPrecision(2)}, ${viewDirection.y.toPrecision(2)}, ${viewDirection.z.toPrecision(2)})`);

	if((viewDirection.x == 0.0) && (viewDirection.y == 0.0)) {
		// viewDirection is along z direction
		apertureBasisVector1.crossVectors(viewDirection, new THREE.Vector3(1, 0, 0)).normalize();
	} else {
		// viewDirection is not along z direction
		apertureBasisVector1.crossVectors(viewDirection, new THREE.Vector3(0, 0, 1)).normalize();
	}
	apertureBasisVector2.crossVectors(viewDirection, apertureBasisVector1).normalize();

	raytracingSphereShaderMaterial.uniforms.noOfRays.value = noOfRays;
	raytracingSphereShaderMaterial.uniforms.apertureXHat.value.copy(apertureBasisVector1);
	raytracingSphereShaderMaterial.uniforms.apertureYHat.value.copy(apertureBasisVector2);
	raytracingSphereShaderMaterial.uniforms.apertureRadius.value = apertureRadius;
	raytracingSphereShaderMaterial.uniforms.focusDistance.value = focusDistance;
}

/** create raytracing phere */
function addRaytracingSphere() {
	const videoFeedUTexture = new THREE.VideoTexture( videoFeedU );
	const videoFeedETexture = new THREE.VideoTexture( videoFeedE );
	videoFeedUTexture.colorSpace = THREE.SRGBColorSpace;
	videoFeedETexture.colorSpace = THREE.SRGBColorSpace;

	// create arrays of random numbers (as GLSL is rubbish at doing random numbers)
	let randomNumbersX = [];
	let randomNumbersY = [];
	// make the first random number 0 in both arrays, meaning the 0th ray starts from the centre of the aperture
	randomNumbersX.push(0);
	randomNumbersY.push(0);
	// fill in the rest of the array with random numbers
	let i=1;
	do {
		// create a new pairs or random numbers (x, y) such that x^2 + y^2 <= 1
		let x = 2*Math.random()-1;	// random number between -1 and 1
		let y = 2*Math.random()-1;	// random number between -1 and 1
		if(x*x + y*y <= 1) {
			// (x,y) lies within a circle of radius 1
			//  add a new point to the array of points on the aperture
			randomNumbersX.push(x);
			randomNumbersY.push(y);
			i++;
		}
	} while (i < 100);

	// the sphere surrouning the camera in all directions
	const geometry = 
		new THREE.SphereGeometry( raytracingSphereRadius );
	raytracingSphereShaderMaterial = new THREE.ShaderMaterial({
		side: THREE.DoubleSide,
		// wireframe: true,
		uniforms: { 
			period: { value: 0.4 },
			focussingPowerGradient: { value: 1 },
			alvarezFocusing: { value: true },
			radius: { value: 5.0 },	// radius of arrays
			visible1: { value: true },
			cosAlpha1: { value: 1.0 },
			sinAlpha1: { value: 0.0	},
			additionalF1: { value: 1e10 },	// additional focal length of lenslet array (an additional lens in the same plane)
			centreOfArray1: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			visible2: { value: true },
			cosAlpha2: { value: 1.0 },
			sinAlpha2: { value: 0.0	},
			additionalF2: { value: 1e10 },	// additional focal length of lenslet array (an additional lens in the same plane)
			centreOfArray2: { value: new THREE.Vector3(0, 0, 0) },	// principal point of lenslet (0, 0)
			videoFeedUTexture: { value: videoFeedUTexture }, 
			videoFeedETexture: { value: videoFeedETexture },
			tanHalfFovHU: { value: 1.0 },
			tanHalfFovVU: { value: 1.0 },
			tanHalfFovHE: { value: 1.0 },
			tanHalfFovVE: { value: 1.0 },
			halfWidthU: { value: 1.0 },
			halfHeightU: { value: 1.0 },
			halfWidthE: { value: 1.0 },
			halfHeightE: { value: 1.0 },
			videoDistance: { value: 10.0 },	// distance of the image of the video feed from the origin
			focusDistance: { value: 10.0 },
			apertureXHat: { value: new THREE.Vector3(1, 0, 0) },
			apertureYHat: { value: new THREE.Vector3(0, 1, 0) },
			apertureRadius: { value: apertureRadius },
			randomNumbersX: { value: randomNumbersX },
			randomNumbersY: { value: randomNumbersY },
			noOfRays: { value: 1 }
		},
		vertexShader: `
			varying vec3 intersectionPoint;

			void main()	{
				// projectionMatrix, modelViewMatrix, position -> passed in from Three.js
				intersectionPoint = position.xyz;
  				gl_Position = projectionMatrix
					* modelViewMatrix
					* vec4(position, 1.0);
			}
		`,
		fragmentShader: `
			precision highp float;

			varying vec3 intersectionPoint;
			
			uniform float period;	// period of both arrays
			uniform float radius;	// radius of both arrays
			uniform float focussingPowerGradient;	// focussing-power gradient dP_x/dx = dP_y/dyd
			uniform bool alvarezFocusing;

			// lenslet array 1
			uniform bool visible1;
			uniform float cosAlpha1;	// cos of rotation angle of array 1
			uniform float sinAlpha1;	// sin of rotation angle of array 1
			uniform float additionalF1;	// additional focal length of lenslet array (an additional lens in the same plane)
			uniform vec3 centreOfArray1;	// centre of array 1,  and principal point of lenslet (0, 0)

			// lenslet array 2
			uniform bool visible2;
			uniform float cosAlpha2;	// cos of rotation angle of array 2
			uniform float sinAlpha2;	// sin of rotation angle of array 2
			uniform float additionalF2;	// additional focal length of lenslet array (an additional lens in the same plane)
			uniform vec3 centreOfArray2;	// centre of array 2,  and principal point of lenslet (0, 0)
			uniform float radius2;	// radius of array 2

			uniform bool idealLenses;	// true => use ideal thin lenses; false => use lens holograms

			// video feed from user-facing camera
			uniform sampler2D videoFeedUTexture;
			uniform float halfWidthU;
			uniform float halfHeightU;

			// video feed from environment-facing camera
			uniform sampler2D videoFeedETexture;
			uniform float halfWidthE;
			uniform float halfHeightE;

			// the camera's wide aperture
			uniform float videoDistance;
			uniform float focusDistance;
			uniform int noOfRays;
			uniform vec3 apertureXHat;
			uniform vec3 apertureYHat;
			uniform float apertureRadius;
			uniform float randomNumbersX[100];
			uniform float randomNumbersY[100];
			// uniform float apertureRadius;

			// rotate the 2D vector v by the angle alpha (in radians)
			// from https://gist.github.com/yiwenl/3f804e80d0930e34a0b33359259b556c
			vec2 rotate(vec2 v, float cosAlpha, float sinAlpha) {
				mat2 m = mat2(cosAlpha, sinAlpha, -sinAlpha, cosAlpha);
				return m * v;
			}

			// propagate the ray starting at position p and with direction d to the plane z = z0, providing that plane
			// is in the ray's "forward" direction;
			// p becomes the point where the ray intersects p;
			// isForward is set to true or false depending on whether the intersection point is forwards or backwards along the ray
			void propagateForwardToZPlane(
				inout vec3 p, 
				vec3 d, 
				float z0,
				inout bool isForward
			) {
				// calculate the z distance from the ray start position to the array
				float deltaZ = z0 - p.z;

				// is the intersection with the plane in the ray's "forward" direction?
				isForward = (d.z*deltaZ > 0.0);

				// if the intersection is in the forward direction, advance the ray to the plane
				if(isForward) p += d/d.z*deltaZ;	// set p to the intersection point with the plane
			}

			// Calculate the light-ray direction after transmission through a lens or lens hologram.
			// d is the incident light-ray direction;
			// pixy is a 2D vector containing the transverse (x, y) components of the vector I-P,
			// i.e. the vector from the principal point P to the intersection point I;
			// f is the focal length;
			// returns the outgoing light-ray direction
			vec3 lensDeflect(vec3 d, vec2 pixy, float f) {
				if(idealLenses) {
					// ideal thin lens
					// "normalise" the direction such that the magnitude of the z component is 1
					vec3 d1 = d/abs(d.z);

					// the 3D deflected direction comprises the transverse components and a z component of magnitude 1
					// and the same sign as d.z
					return vec3(d1.xy - pixy/f, d1.z);
				} else {
					// lens hologram
					// normalise d
					vec3 dN = d/length(d);
					// transverse components of the outgoing light-ray direction
					vec2 dxy = dN.xy - pixy/f;
	
					// from the transverse direction, construct a 3D vector by setting the z component such that the length
					// of the vector is 1
					return vec3(dxy, sign(d.z)*sqrt(1.0 - dot(dxy, dxy)));
				}
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) a lens.
			// The (ideal thin) lens, of focal length f, is in a z plane through centreOfLens.
			// It is circular, with the given radius, centred on centreOfLenss.
			void passThroughLens(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				vec3 centreOfLens, 
				float radius,
				float focalLength
			) {
				bool isForward;
				propagateForwardToZPlane(p, d, centreOfLens.z, isForward);

				if(isForward) {
					// there is an intersection with the plane of this lens in the ray's forward direction

					// does the intersection point lie within the radius?
					vec2 pixy = p.xy - centreOfLens.xy;
					float r2 = dot(pixy, pixy);
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the lens does something to the ray

						// deflect the light-ray direction accordingly and make sure that the sign of the z component remains the same
						lensDeflect(d, pixy, focalLength);

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} 
				}
			}

			float findLensletCentreCoordinate(float u, float uPeriod) {
				return uPeriod*floor(u/uPeriod+0.5);
			}

			vec3 crossedLinearPowerLenticularArraysDeflect(vec3 d, vec3 intersectionPoint, vec3 principalPoint00, float cosAlpha, float sinAlpha, float dPdUV) {
				vec2 r = intersectionPoint.xy - principalPoint00.xy;

				// transform r into the uv coordinate system
				vec2 rUV = rotate(r, cosAlpha, -sinAlpha);

				// lenslet centre, in the uv coordinate system aligned with the array vectors
				vec2 lensletCentreUV = vec2(
					findLensletCentreCoordinate(rUV.x, period),
					findLensletCentreCoordinate(rUV.y, period)
				);

				// calculate the tangential direction change
				vec2 localRUV = rUV - lensletCentreUV;
				vec2 deltaUV = -localRUV*dPdUV*rUV;
				if(alvarezFocusing) deltaUV += 0.5*dPdUV*localRUV*localRUV;	// why "+"?
				// the second term, proportional to -0.5*dPdUV*uuN^2 (or vvN^2) is for one part of a Lohmann cylindrical lens --
				// see PhaseHologramOfLohmannCylindricalLensPart::getTangentialDirectionComponentChangeTransmissive
		
				// transform the lenslet centre into the xy coordinate system
				vec2 dxy = d.xy/length(d) + rotate(deltaUV, cosAlpha, sinAlpha);

				// from the transverse direction, construct a 3D vector by setting the z component such that the length
				// of the vector is 1
				return vec3(dxy, sign(d.z)*sqrt(1.0 - dot(dxy, dxy)));
			}

			// Pass the current ray (start point p, direction d, brightness factor b) through (or around) crossed
			// linear-power lenticular arrays.
			// The lenslet array is in a z plane through centreOfArray; it is simulated as phase holograms.
			// The component is essentially the sum of two crossed lenticular arrays of the same period.
			// The focal power in both arrays is zero in the centre and it is proportional to position,
			// such that P_u = dPdUV * u and P_v = dPdUV * v.
			// The whole array isalso rotated by alpha around the z axis.  
			// The component is circular, with the given radius, centred on centreOfArray.
			void passThroughCrossedLinearPowerLenticularArrays(
				inout vec3 p, 
				inout vec3 d, 
				inout vec4 b,
				vec3 centreOfArray, 
				float radius,
				float cosAlpha,
				float sinAlpha,
				float dPdUV,
				float overallF
			) {
				bool isForward;
				propagateForwardToZPlane(p, d, centreOfArray.z, isForward);

				if(isForward) {
					// there is an intersection with the plane of this array in the ray's forward direction
					
					// does the intersection point lie with in the radius?
					vec2 pixy = p.xy - centreOfArray.xy;
					float r2 = dot(pixy, pixy);	// length squared of vector r
					if(r2 < radius*radius) {
						// the intersection point lies inside the radius, so the component does something to the ray

						// deflect the light-ray direction accordingly 

						d = crossedLinearPowerLenticularArraysDeflect(d, p, centreOfArray, cosAlpha, sinAlpha, dPdUV);
						d = lensDeflect(d, pixy, overallF);

						// lower the brightness factor, giving the light a blue tinge
						b *= vec4(0.9, 0.9, 0.99, 1);
					} 
				} else b *= vec4(0.99, 0.9, 0.9, 1);	// this shouldn't happen -- give the light a red tinge
			}

			// propagate the ray to the plane of the video feed, which is a z-distance <videoDistance> away,
			// and return either the color of the corresponding video-feed texel or the background color
			vec4 getColorOfVideoFeed(
				inout vec3 p, 
				vec3 d, 
				vec4 b,
				float videoFeedZ,
				sampler2D videoFeedTexture,
				float halfWidth,
				float halfHeight,
				vec4 backgroundColor
			) {
				bool isForward;
				propagateForwardToZPlane(p, d, videoFeedZ, isForward);

				// is the intersection in the ray's forward direction?
				if(isForward) {
					// does the ray intersect the image?
					if((abs(p.x) < halfWidth) && (abs(p.y) < halfHeight))
						// yes, the ray intersects the image; take the pixel colour from the camera's video feed
						return texture2D(videoFeedTexture, vec2(0.5+0.5*p.x/halfWidth, 0.5+0.5*p.y/halfHeight));
					else 
						// the ray doesn't intersect the image
						return backgroundColor;
				}
			}

			void main() {
				// first calculate the point this pixel is focussed on, which is in a z plane a distance
				// <focusDistance> in  front of the camera, in the direction from the camera's aperture centre to the intersection point
				vec3 dF = intersectionPoint - cameraPosition;
				vec3 focusPosition = cameraPosition + focusDistance/abs(dF.z)*dF;

				// trace <noOfRays> rays
				gl_FragColor = vec4(0, 0, 0, 0);
				vec4 color;
				for(int i=0; i<noOfRays; i++) {
					// the current ray start position, a random point on the camera's circular aperture
					vec3 p = cameraPosition + apertureRadius*randomNumbersX[i]*apertureXHat + apertureRadius*randomNumbersY[i]*apertureYHat;
	
					// first calculate the current light-ray direction:
					// the ray first passes through focusPosition and then p,
					// so the "backwards" ray direction from the camera to the intersection point is
					//   d = focusPosition - p
					vec3 d = focusPosition - p;
					d = dF.z/d.z*d;
	
					// current brightness factor; this will multiply the colour at the end
					vec4 b = vec4(1.0, 1.0, 1.0, 1.0);
	
					if(d.z < 0.0) {
						// the ray is travelling "forwards", in the (-z) direction;
						// pass first through array 1, then array 2, then to environment-facing video feed
						if(visible1) passThroughCrossedLinearPowerLenticularArrays(p, d, b, centreOfArray1, radius, cosAlpha1, sinAlpha1, +focussingPowerGradient, additionalF1);
						if(visible2) passThroughCrossedLinearPowerLenticularArrays(p, d, b, centreOfArray2, radius, cosAlpha2, sinAlpha2, -focussingPowerGradient, additionalF2);
						color = getColorOfVideoFeed(p, d, b, -videoDistance, videoFeedETexture, halfWidthE, halfHeightE, vec4(1, 1, 1, 1.0));	// white
					} else {
						// the ray is travelling "backwards", in the (+z) direction;
						// pass first through array 2, then array 1, then to user-facing video feed
						if(visible2) passThroughCrossedLinearPowerLenticularArrays(p, d, b, centreOfArray2, radius, cosAlpha2, sinAlpha2, -focussingPowerGradient, additionalF2);
						if(visible1) passThroughCrossedLinearPowerLenticularArrays(p, d, b, centreOfArray1, radius, cosAlpha1, sinAlpha1, +focussingPowerGradient, additionalF1);
						color = getColorOfVideoFeed(p, d, b, videoDistance, videoFeedUTexture, halfWidthU, halfHeightU, vec4(1, 0, 0, 1.0));	// white
					}
		
					// finally, multiply by the brightness factor and add to gl_FragColor
					gl_FragColor += b*color;
				}
					
				gl_FragColor /= float(noOfRays);
			}
		`
	});
	raytracingSphere = new THREE.Mesh( geometry, raytracingSphereShaderMaterial ); 
	scene.add( raytracingSphere );
}

function createVideoFeeds() {
	// create the video stream for the user-facing camera first, as some devices (such as my iPad), which have both cameras,
	// but can (for whatever reason) only have a video feed from one at a time, seem to go with the video stream that was
	// created last, and as the standard view is looking "forward" it is preferable to see the environment-facing camera.
	videoFeedU = document.getElementById( 'videoFeedU' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// user-facing camera
		const constraintsU = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'user'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraintsU ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeedU.srcObject = stream;
			videoFeedU.play();

			videoFeedU.addEventListener("playing", () => {
				aspectRatioVideoFeedU = videoFeedU.videoWidth / videoFeedU.videoHeight;
				updateUniforms();
				postStatus(`User-facing(?) camera resolution ${videoFeedU.videoWidth} &times; ${videoFeedU.videoHeight}`);
			});
		} ).catch( function ( error ) {
			postStatus(`Unable to access user-facing camera/webcam (Error: ${error})`);
		} );
	} else {
		postStatus( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}

	videoFeedE = document.getElementById( 'videoFeedE' );

	// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_materials_video_webcam.html
	if ( navigator.mediaDevices && navigator.mediaDevices.getUserMedia ) {
		// environment-facing camera
		const constraintsE = { video: { 
			// 'deviceId': cameraId,	// this could be the device ID selected 
			width: {ideal: 1280},	// {ideal: 10000}, 
			// height: {ideal: 10000}, 
			facingMode: {ideal: 'environment'}
			// aspectRatio: { exact: width / height }
		} };
		navigator.mediaDevices.getUserMedia( constraintsE ).then( function ( stream ) {
			// apply the stream to the video element used in the texture
			videoFeedE.srcObject = stream;
			videoFeedE.play();

			videoFeedE.addEventListener("playing", () => {
				aspectRatioVideoFeedE = videoFeedE.videoWidth / videoFeedE.videoHeight;
				updateUniforms();
				postStatus(`Environment-facing(?) camera resolution ${videoFeedE.videoWidth} &times; ${videoFeedE.videoHeight}`);
			});
		} ).catch( function ( error ) {
			postStatus(`Unable to access environment-facing camera/webcam (Error: ${error})`);
		} );
	} else {
		postStatus( 'MediaDevices interface, which is required for video streams from device cameras, not available.' );
	}
}

function addEventListenersEtc() {
	// handle device orientation
	// window.addEventListener("deviceorientation", handleOrientation, true);

	// handle window resize
	window.addEventListener("resize", onWindowResize, false);

	// handle screen-orientation (landscape/portrait) change
	screen.orientation.addEventListener( "change", recreateVideoFeeds );

	// share button functionality
	document.getElementById('takePhotoButton').addEventListener('click', takePhoto);

	// toggle fullscreen button functionality
	document.getElementById('fullscreenButton').addEventListener('click', toggleFullscreen);

	// info button functionality
	document.getElementById('infoButton').addEventListener('click', toggleInfoVisibility);

	// back button functionality
	document.getElementById('backButton').addEventListener('click', showLivePhoto);
	document.getElementById('backButton').style.visibility = "hidden";

	// share button
	document.getElementById('shareButton').addEventListener('click', share);
	document.getElementById('shareButton').style.visibility = "hidden";
	if(!(navigator.share)) document.getElementById('shareButton').src="./shareButtonUnavailable.png";
	// if(!(navigator.share)) document.getElementById('shareButton').style.opacity = 0.3;

	// delete button
	document.getElementById('deleteButton').addEventListener('click', deleteStoredPhoto);
	document.getElementById('deleteButton').style.visibility = "hidden";

	// hide the thumbnail for the moment
	document.getElementById('storedPhotoThumbnail').addEventListener('click', showStoredPhoto);
	document.getElementById('storedPhotoThumbnail').style.visibility = "hidden";
	document.getElementById('storedPhoto').addEventListener('click', showLivePhoto);
	document.getElementById('storedPhoto').style.visibility = "hidden";
	// showingStoredPhoto = false;
}

// see https://github.com/mrdoob/three.js/blob/master/examples/webgl_animation_skinning_additive_blending.html
function createGUI() {
	// const 
	gui = new GUI();
	// gui.hide();

	const params = {
		'&Delta;<i>x</i>': offset.x,
		'&Delta;<i>y</i>': offset.y,
		'Show comp. 1': raytracingSphereShaderMaterial.uniforms.visible1.value,
		'Show comp. 2': raytracingSphereShaderMaterial.uniforms.visible2.value,
		'd<i>P<sub>x</sub></i>/d<i>x</i> = d<i>P<sub>y</sub></i>/d<i>y</i>': raytracingSphereShaderMaterial.uniforms.focussingPowerGradient.value,
		'Period': raytracingSphereShaderMaterial.uniforms.period.value,
		'Alvarez focussing': raytracingSphereShaderMaterial.uniforms.alvarezFocusing.value,
		'Radius': raytracingSphereShaderMaterial.uniforms.radius.value,	// radius of the Fresnel lens
		'&Delta;<i>z</i>': offset.z,
		'&Delta;&phi; (&deg;)': deltaPhi / Math.PI * 180.,
		'1/<i>f</i><sub>1</sub>': 0.0,	
		'1/<i>f</i><sub>2</sub>': 0.0,	
		// 'Show equivalent ideal lens': raytracingSphereShaderMaterial.uniforms.showEquivalentLens.value,
		'Horiz. FOV (&deg;)': fovScreen,
		'Aperture radius': apertureRadius,
		'1/(focus. dist.)': 1.0/focusDistance,
		'No of rays': noOfRays,
		'Env.-facing cam. (&deg;)': fovVideoFeedE,
		'User-facing cam. (&deg;)': fovVideoFeedU,
		'tan<sup>-1</sup>(video dist.)': Math.atan(raytracingSphereShaderMaterial.uniforms.videoDistance.value),
		'Point (virtual) cam. forward (in -<b>z</b> direction)': pointForward,
		'Show/hide info': toggleInfoVisibility,
		'Restart video streams': function() { 
			recreateVideoFeeds(); 
			postStatus("Restarting video stream");
		}
	}

	gui.add( params, '&Delta;<i>x</i>', -0.1, 0.1).onChange( (d) => { offset.x = d; } );
	gui.add( params, '&Delta;<i>y</i>', -0.1, 0.1).onChange( (d) => { offset.y = d; } );

	const otherArrayParams = gui.addFolder( 'Components' );

	otherArrayParams.add( params, 'Show comp. 1').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible1.value = v; } );
	otherArrayParams.add( params, 'Show comp. 2').onChange( (v) => { raytracingSphereShaderMaterial.uniforms.visible2.value = v; } );
	otherArrayParams.add( params, 'd<i>P<sub>x</sub></i>/d<i>x</i> = d<i>P<sub>y</sub></i>/d<i>y</i>', -10, 10).onChange( (dPxydxy) => { raytracingSphereShaderMaterial.uniforms.focussingPowerGradient.value = dPxydxy; } );
	otherArrayParams.add( params, 'Alvarez focussing' ).onChange( (a) => { raytracingSphereShaderMaterial.uniforms.alvarezFocusing.value = a; } );
	otherArrayParams.add( params, 'Period', 0.01, 1).onChange( (p) => { raytracingSphereShaderMaterial.uniforms.period.value = p; } );
	otherArrayParams.add( params, '1/<i>f</i><sub>1</sub>', -1, 1).onChange( (p) => { raytracingSphereShaderMaterial.uniforms.additionalF1.value = 1.0/p; } );
	otherArrayParams.add( params, '1/<i>f</i><sub>2</sub>', -1, 1).onChange( (p) => { raytracingSphereShaderMaterial.uniforms.additionalF2.value = 1.0/p; } );
	otherArrayParams.add( params, 'Radius', 0.1, 10 ).onChange( (r) => {raytracingSphereShaderMaterial.uniforms.radius.value = r; } );
	
	otherArrayParams.add( params, '&Delta;<i>z</i>', 0.00001, 1).onChange( (d) => { offset.z = d; } );
	otherArrayParams.add( params, '&Delta;&phi; (&deg;)', -10, 10).onChange( (alpha) => { deltaPhi = alpha/180.0*Math.PI; } );
	// folderArray1.add( params, 'Show equivalent ideal lens' ).onChange( (s) => {raytracingSphereShaderMaterial.uniforms.showEquivalentLens.value = s; } );

	const folderDevice = gui.addFolder( 'Device cameras horiz. FOV' );
	folderDevice.add( params, 'Env.-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedE = fov; updateUniforms(); });   
	folderDevice.add( params, 'User-facing cam. (&deg;)', 10, 170, 1).onChange( (fov) => { fovVideoFeedU = fov; updateUniforms(); });   
	folderDevice.close();

	const folderVirtualCamera = gui.addFolder( 'Virtual camera' );
	folderVirtualCamera.add( params, 'Horiz. FOV (&deg;)', 10, 170, 1).onChange( setScreenFOV );
	folderVirtualCamera.add( params, 'Aperture radius', 0.0, 1.0).onChange( (r) => { apertureRadius = r; } );
	folderVirtualCamera.add( params, '1/(focus. dist.)', -1, 1).onChange( (i) => { focusDistance = 1.0/i; } );
	folderVirtualCamera.add( params, 'No of rays', 1, 100, 1).onChange( (n) => { noOfRays = n; } );

	const folderSettings = gui.addFolder( 'Other controls' );
	folderSettings.add( params, 'tan<sup>-1</sup>(video dist.)', Math.atan(0.1), 0.5*Math.PI).onChange( (a) => { raytracingSphereShaderMaterial.uniforms.videoDistance.value = Math.tan(a); } );
	folderSettings.add( params, 'Lenslet type', { 'Ideal thin': true, 'Phase hologram': false } ).onChange( (t) => { raytracingSphereShaderMaterial.uniforms.idealLenses.value = t; });
	// folderSettings.add( params, 'Ideal lenses').onChange( (b) => { raytracingSphereShaderMaterial.uniforms.idealLenses.value = b; } );
	folderSettings.add( params, 'Point (virtual) cam. forward (in -<b>z</b> direction)');
	folderSettings.add( params, 'Show/hide info');
	folderSettings.add( params, 'Restart video streams');
	folderSettings.close();
}

/**
 * @param {*} fov	The larger of the camera's horizontal and vertical FOV, in degrees
 * 
 * Set the larger FOV of the screen/window to fov.
 * 
 * Depending on the screen/window's FOV, fov is either the horizontal fov (if screen width > screen height)
 * or the vertical fov (if screen width < screen height).
 */
function setScreenFOV(fov) {
	fovScreen = fov;

	screenChanged();
}

// function swapArrays() {
// 	const visible3 = lensletArrayShaderMaterial.uniforms.visible1.value;
// 	const focalLength3 = lensletArrayShaderMaterial.uniforms.lensletsF1.value;
// 	const period3 = lensletArrayShaderMaterial.uniforms.period1.value;
// 	const alpha3 = lensletArrayShaderMaterial.uniforms.alpha1.value;

// 	lensletArrayShaderMaterial.uniforms.visible1.value = lensletArrayShaderMaterial.uniforms.visible2.value;
// 	lensletArrayShaderMaterial.uniforms.lensletsF1.value = lensletArrayShaderMaterial.uniforms.lensletsF2.value;
// 	lensletArrayShaderMaterial.uniforms.period1.value = lensletArrayShaderMaterial.uniforms.period2.value;
// 	lensletArrayShaderMaterial.uniforms.alpha1.value = lensletArrayShaderMaterial.uniforms.alpha2.value;

// 	lensletArrayShaderMaterial.uniforms.visible2.value = visible3;
// 	lensletArrayShaderMaterial.uniforms.lensletsF2.value = focalLength3;
// 	lensletArrayShaderMaterial.uniforms.period2.value = period3;
// 	lensletArrayShaderMaterial.uniforms.alpha2.value = alpha3;
// }

/** 
 * Reset the aspect ratio and FOV of the virtual cameras.
 * 
 * Call if the window size has changed (which also happens when the screen orientation changes)
 * or if camera's FOV has changed
 */
function screenChanged() {
	// alert(`new window size ${window.innerWidth} x ${window.innerHeight}`);

	// in case the screen size has changed
	if(renderer) renderer.setSize(window.innerWidth, window.innerHeight);

	// if the screen orientation changes, width and height swap places, so the aspect ratio changes
	let windowAspectRatio = window.innerWidth / window.innerHeight;
	camera.aspect = windowAspectRatio;

	// fovS is the screen's horizontal or vertical FOV, whichever is greater;
	// re-calculate the camera FOV, which is the *vertical* fov
	let verticalFOV;
	if(windowAspectRatio > 1.0) {
		// fovS is horizontal FOV; convert to get correct vertical FOV
		verticalFOV = 2.0*Math.atan(Math.tan(0.5*fovScreen*Math.PI/180.0)/windowAspectRatio)*180.0/Math.PI;
	} else {
		// fovS is already vertical FOV
		verticalFOV = fovScreen;
	}
	camera.fov = verticalFOV;

	// make sure the camera changes take effect
	camera.updateProjectionMatrix();
}

function  pointForward() {
	let r = camera.position.length();
	camera.position.x = 0;
	camera.position.y = 0;
	camera.position.z = r;
	controls.update();
	postStatus('Pointing camera forwards (in -<b>z</b> direction)');
}

function onWindowResize() {
	screenChanged();
	postStatus(`window size ${window.innerWidth} &times; ${window.innerHeight}`);	// debug
}

// // see https://developer.mozilla.org/en-US/docs/Web/API/ScreenOrientation/change_event
function recreateVideoFeeds() {
	// stop current video streams...
	videoFeedE.srcObject.getTracks().forEach(function(track) { track.stop(); });
	videoFeedU.srcObject.getTracks().forEach(function(track) { track.stop(); });

	// ... and re-create new ones, hopefully of the appropriate size
	createVideoFeeds();
}

function addOrbitControls() {
	// controls

	controls = new OrbitControls( camera, renderer.domElement );
	// controls = new OrbitControls( cameraOutside, renderer.domElement );
	controls.listenToKeyEvents( window ); // optional

	//controls.addEventListener( 'change', render ); // call this only in static scenes (i.e., if there is no animation loop)
	controls.addEventListener( 'change', cameraPositionChanged );

	controls.enableDamping = false; // an animation loop is required when either damping or auto-rotation are enabled
	controls.dampingFactor = 0.05;

	controls.enablePan = true;
	controls.enableZoom = true;

	controls.maxPolarAngle = Math.PI;
}

function cameraPositionChanged() {
	postStatus(`Camera position (${camera.position.x.toPrecision(2)}, ${camera.position.y.toPrecision(2)}, ${camera.position.z.toPrecision(2)})`);
	// counter = 0;
	// keep the raytracing sphere centred on the camera position
	// raytracingSphere.position.copy(camera.position.clone());	// TODO this doesn't seem to work as intended!?
}

async function toggleFullscreen() {
	if (!document.fullscreenElement) {
		document.documentElement.requestFullscreen().catch((err) => {
			postStatus(
				`Error attempting to enable fullscreen mode: ${err.message} (${err.name})`,
			);
		});
		// allow screen orientation changes
		// screen.orientation.unlock();
	} else {
		document.exitFullscreen();
	}
}

function showStoredPhoto() {
	gui.hide();
	renderer.domElement.style.visibility = "hidden";
	document.getElementById('takePhotoButton').style.visibility = "hidden";
	// document.getElementById('changePositionButton').style.visibility = "hidden";
	document.getElementById('storedPhotoThumbnail').style.visibility = "hidden";
	document.getElementById('backButton').style.visibility = "visible";
	document.getElementById('shareButton').style.visibility = "visible";
	document.getElementById('deleteButton').style.visibility = "visible";
	document.getElementById('storedPhoto').style.visibility = "visible";
	showingStoredPhoto = true;

	postStatus('Showing stored photo, '+storedPhotoDescription);
}

function showLivePhoto() {
	gui.show();
	renderer.domElement.style.visibility = "visible";
	document.getElementById('takePhotoButton').style.visibility = "visible";
	// document.getElementById('changePositionButton').style.visibility = "visible";
	if(storedPhoto) document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	document.getElementById('backButton').style.visibility = "hidden";
	document.getElementById('shareButton').style.visibility = "hidden";
	document.getElementById('deleteButton').style.visibility = "hidden";
	document.getElementById('storedPhoto').style.visibility = "hidden";
	showingStoredPhoto = false;

	postStatus('Showing live image');
}

function deleteStoredPhoto() {
	storedPhoto = null;

	showLivePhoto();

	postStatus('Stored photo deleted; showing live image');
}

function takePhoto() {
	try {
		click.play();

		storedPhoto = renderer.domElement.toDataURL('image/png');
		storedPhotoInfoString = getInfoString();

		storedPhotoDescription = name;
		// 
		document.getElementById('storedPhoto').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').src=storedPhoto;
		document.getElementById('storedPhotoThumbnail').style.visibility = "visible";
	
		postStatus('Photo taken; click thumbnail to view and share');
	} catch (error) {
		console.error('Error:', error);
	}	
}

async function share() {
	try {
		fetch(storedPhoto)
		.then(response => response.blob())
		.then(blob => {
			const file = new File([blob], name+storedPhotoDescription+'.png', { type: blob.type });

			// Use the Web Share API to share the screenshot
			if (navigator.share) {
				navigator.share({
					title: storedPhotoDescription,
					files: [file],
				});
			} else {
				postStatus('Sharing is not supported by this browser.');
			}	
		})
		.catch(error => {
			console.error('Error:', error);
			postStatus(`Error: ${error}`);
		});
	} catch (error) {
		console.error('Error:', error);
	}
}

/** 
 * Add a text field to the bottom left corner of the screen
 */
function createStatus() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	status.style.position = 'absolute';
	status.style.backgroundColor = "rgba(0, 0, 0, 0.3)";	// semi-transparent black
	status.style.color = "White";
	status.style.fontFamily = "Arial";
	status.style.fontSize = "9pt";
	postStatus("Welcome!");
	status.style.bottom = 0 + 'px';
	status.style.left = 0 + 'px';
	status.style.zIndex = 1;
	document.body.appendChild(status);	
}

function postStatus(text) {
	status.innerHTML = '&nbsp;'+text;
	console.log('status: '+text);

	// show the text only for 3 seconds
	statusTime = new Date().getTime();
	setTimeout( () => { if(new Date().getTime() - statusTime > 2999) status.innerHTML = '&nbsp;'+name+', University of Glasgow, <a href="https://github.com/jkcuk/'+name+'">https://github.com/jkcuk/'+name+'</a>' }, 3000);
}

function getInfoString() {
	return `Lenslet arrays<br>` +
		`&nbsp;&nbsp;Period = ${raytracingSphereShaderMaterial.uniforms.period.value.toPrecision(4)}<br>` +
		'&nbsp;&nbsp;Alvarez focussing ' + (raytracingSphereShaderMaterial.uniforms.alvarezFocusing.value?'&check;':'&cross;')+`<br>\n` +
		`&nbsp;&nbsp;Clear-aperture radius = ${raytracingSphereShaderMaterial.uniforms.radius.value.toPrecision(4)}<br>\n` +	// radius of the Fresnel lens
		`&nbsp;&nbsp;Rotation angle &Delta;&phi; = ${(deltaPhi*180.0/Math.PI).toPrecision(4)}&deg;<br>` +
		`&nbsp;&nbsp;Lenslet array 1 (the closer array, when seen in "forward" direction)<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;Visible `+ (raytracingSphereShaderMaterial.uniforms.visible1.value?'&check;':'&cross;')+`<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;Centre of array = (${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray1.value.z.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;Focal length of additional lens in same plane = ${raytracingSphereShaderMaterial.uniforms.additionalF1.value.toPrecision(4)}<br>` +		
		`&nbsp;&nbsp;Lenslet array 2 (the farther array, when seen in "forward" direction)<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;Visible `+ (raytracingSphereShaderMaterial.uniforms.visible2.value?'&check;':'&cross;')+`<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;Centre of array = (${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.centreOfArray2.value.z.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;Focal length of additional lens in same plane = ${raytracingSphereShaderMaterial.uniforms.additionalF2.value.toPrecision(4)}<br>` +		
		`Video feeds<br>` +
		`&nbsp;&nbsp;Distance from origin = ${raytracingSphereShaderMaterial.uniforms.videoDistance.value.toPrecision(4)}<br>` +	// (user-facing) camera
		`&nbsp;&nbsp;Horizontal fields of view (when seen from the origin)<br>` +
		`&nbsp;&nbsp;&nbsp;&nbsp;User-facing camera = ${fovVideoFeedU.toPrecision(4)}&deg;<br>` +	// (user-facing) camera
		`&nbsp;&nbsp;&nbsp;&nbsp;Environment-facing camera = ${fovVideoFeedE.toPrecision(4)}&deg;<br>` +	// (environment-facing) camera
		`Virtual camera<br>` +
		`&nbsp;&nbsp;Position = (${camera.position.x.toPrecision(4)}, ${camera.position.y.toPrecision(4)}, ${camera.position.z.toPrecision(4)})<br>` +
		`&nbsp;&nbsp;Horiz. FOV = ${fovScreen.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Aperture radius = ${apertureRadius.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Focussing distance = ${focusDistance.toPrecision(4)}<br>` +
		`&nbsp;&nbsp;Number of rays = ${noOfRays}`
		// `apertureXHat = (${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureXHat.value.z.toPrecision(4)})<br>` +
		// `apertureYHat = (${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.x.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.y.toPrecision(4)}, ${raytracingSphereShaderMaterial.uniforms.apertureYHat.value.z.toPrecision(4)})`
		;
		console.log("*");
}

function refreshInfo() {
	if(showingStoredPhoto) setInfo( storedPhotoInfoString );
	else setInfo( getInfoString() );

	if(info.style.visibility == "visible") setTimeout( refreshInfo , 100);	// refresh again a while
}

/** 
 * Add a text field to the top left corner of the screen
 */
function createInfo() {
	// see https://stackoverflow.com/questions/15248872/dynamically-create-2d-text-in-three-js
	info.style.position = 'absolute';
	info.style.backgroundColor = "rgba(0, 0, 0, 0.3)";	// semi-transparent black
	info.style.color = "White";
	info.style.fontFamily = "Arial";
	info.style.fontSize = "9pt";
	info.innerHTML = "-- nothing to show (yet) --";
	info.style.top = 60 + 'px';
	info.style.left = 0 + 'px';
	info.style.zIndex = 1;
	document.body.appendChild(info);
	info.style.visibility = "hidden";
}

function setInfo(text) {
	info.innerHTML = text;
	console.log('info: '+text);
	// // show the text only for 3 seconds
	// infoTime = new Date().getTime();
	// setTimeout( () => { if(new Date().getTime() - infoTime > 2999) info.innerHTML = `` }, 3000);
	// info.style.visibility = "visible";
}

function toggleInfoVisibility() {
	switch(info.style.visibility) {
		case "visible":
			info.style.visibility = "hidden";
			break;
		case "hidden":
		default:
			info.style.visibility = "visible";
			refreshInfo();
	}
}