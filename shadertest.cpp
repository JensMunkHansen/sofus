"use strict";

const vs = `#version 300 es

    in vec4 position;

void main() {
  gl_Position = position;
}
`;

const fs = `#version 300 es
    precision mediump float;
precision lowp sampler3D;

uniform sampler3D u_someTexture;

out vec4 theColor;

void main() {
  theColor = texture(u_someTexture, vec3(0,0,0));
}
`;

function main() {
  var m4 = twgl.m4;
  var gl = twgl.getContext(document.createElement("canvas"));
  log("using: " + gl.getParameter(gl.VERSION));
  if (!twgl.isWebGL2(gl)) {
    log("Sorry, this example requires WebGL 2.0");
    return;
  }

  var programInfo = twgl.createProgramInfo(gl, [vs, fs], (err) => {
      log("could not compile shader: " + err);
    });
  if (programInfo) {
    log("compiled shader with sampler3D");
  }

}
main();

function log() {
  var elem = document.createElement("pre");
  elem.appendChild(document.createTextNode(Array.prototype.join.call(arguments, " ")));
  document.body.appendChild(elem);
}
