#!/usr/local/bin/MathematicaScript -runfirst "$TopDirectory=\"/usr/local/Wolfram/Mathematica/10.0\"" -script
symboliclegendre[n_, x_] := Solve[LegendreP[n, x] == 0]; 

legendreprime[n_, a_] := D[LegendreP[n, x], x] /. x -> a; 

weights[n_, x_] := 2/((1 - x^2)*legendreprime[n, x]^2); 

If[Length[$ScriptCommandLine] < 3, Print["Calling:\n",$ScriptCommandLine[[1]]," [OutputName] [Order]"];Exit[]]

h = FromDigits[$ScriptCommandLine[[3]]];
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

WriteString[str,"namespace gl {"];
Write[str];
WriteString[str, "  /* abcissae */"]; 
Do[Write[str]; WriteString[str, "  const static double abcissa", n, "[", Floor[n/2], "]= "]; 
      WriteString[str, "{"]; nlist = symboliclegendre[n, x]; 
  xnlist = x /. nlist; 
      xnlist = N[Abs[Re[xnlist]], precision]; 
  order = Ordering[xnlist]; 
      order = Reverse[order[[1 + Mod[n, 2] ;; n ;; 2]]]; 
  Do[WriteString[str, format[xnlist[[i]]], ","], 
        {i, order[[1 ;; -2]]}]; 
  WriteString[str, format[xnlist[[order[[-1]]]]], "};"]; 
  Write[str]; , {n, 2, h}]; 
Write[str]; 
(* We are ordering the absolute values, which causes an extra reverse *)
WriteString[str, "  /* weights */"]; 
Do[Write[str]; 
  WriteString[str, "  const static double weight", n, "[", Ceiling[n/2], "]= "]; 
      WriteString[str, "{"]; slist := symboliclegendre[n, x]; 
  xslist = x /. slist; 
      xnlist = N[Abs[ArcCos[Abs[Re[xslist]]]], precision]; 
  order = Ordering[xnlist]; 
      order = Reverse[order[[1 ;; n ;; 2]]]; 
      Do[WriteString[str, format[N[Re[weights[n, xslist[[i]]]], precision]], 
    ","], {i, order[[1 ;; -2]]}]; 
      WriteString[str, 
   format[N[Re[weights[n, xslist[[order[[-1]]]]]], precision]], "};"]; 
  Write[str]; , 
     {n, 2, h}]; 

Write[str];
WriteString[str, "  const double* weights[", h - 1, "] = {"]; 
Do[WriteString[str, "weight", n, ", "], {n, 2, h - 1}]; 
WriteString[str, "weight", h, "};"]; 
Write[str];
Write[str];
WriteString[str, "  const double* abcissas[", h - 1, "] = {"]; 
Do[WriteString[str, "abcissa", n, ", "], {n, 2, h - 1}]; 
WriteString[str, "abcissa", h, "};"]; 
Write[str];
WriteString[str,"}"];
Write[str];
Close[str]; 
