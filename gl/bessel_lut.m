#!/usr/local/bin/MathematicaScript -runfirst "$TopDirectory=\"/usr/local/Wolfram/Mathematica/10.0\"" -script

If[Length[$ScriptCommandLine] < 3, Print["Calling:\n",$ScriptCommandLine[[1]]," [OutputName] [Order]"];Exit[]]

n = FromDigits[$ScriptCommandLine[[3]]];
precision = Ceiling[N[Log[10, 2^52]]]; 

format[x_] := PaddedForm[x, {precision, precision-1}, ExponentFunction -> (1 Quotient[#, 1] &), NumberFormat -> (SequenceForm[#1, "e", #3] &)];

str = OpenWrite[$ScriptCommandLine[[2]]];

WriteString[str,"/**\n"];
WriteString[str," * @file   ",$ScriptCommandLine[[2]],"\n"];
WriteString[str," * @author Jens Munk Hansen <jens.munk.hansen@gmai.com>\n"];
WriteString[str," * @date   Tue Apr 29 19:46:27 2014\n"];
WriteString[str," *\n"];
WriteString[str," * @brief  Auto-generated LUT file\n"];
WriteString[str," *\n"];
WriteString[str," *\n"];
WriteString[str," */\n"];
Write[str];
(* k'th zeros of BesselJ(0,z) *)
(* Was const static double JZ when outside namespace *)
WriteString[str,"namespace gl {\n"];
WriteString[str, "  /* k'th zero of BesselJ(0,z) */"]; 
Write[str];
WriteString[str, "  double JZ", "[", n, "]= "]; 
WriteString[str, "{"];
Do[WriteString[str, format[N[BesselJZero[0,i], precision]], ", "]; , {i,1,n-1}];
WriteString[str,format[N[BesselJZero[0,n], precision]], "};\n"];
Write[str];

(* Square of BesselJ(1, BesselZero(0,k)) *)
m = n + 1;
WriteString[str, "  /* Square of BesselJ(1,BesselZero(0,k)) */\n"]; 
WriteString[str, "  double J1", "[", m, "]= "]; 
WriteString[str, "{"];
Do[WriteString[str, format[N[BesselJ[1, BesselJZero[0, i]]^2, precision]], ", "]; , {i,1,m-1}];
WriteString[str, format[N[BesselJ[1, BesselJZero[0, m]]^2, precision]], "};"];
Write[str];
WriteString[str,"}"];
Write[str];
Close[str]; 
