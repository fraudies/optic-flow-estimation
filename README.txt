This directory contains Matlab (R) implementations of methods that estimate 
optic flow from image sequences.

I organized the implementations into folders named according to the authors
of the method. In alphabetical order, these folders are:
AdelsonBergen
Farnebaeck
FleetJepson
Heeger
HornSchunck
LucasKanade
MotaEtAl
Nagel
OtteNagel
ShizawaMase
UrasEtAl

Each folder contains at least the following two scripts (one is a function):
'DemoEstimateOpticFlow2D.m' and 'estiamteOpticFlow2D.m'. The first script
is the method that demonstrates the estimation of optic flow for a particular
method and the second script is the implementation of that method. The demo
script can be run by pressing F5. Note that several methods are compute intensive
and you have to wait 10-60sec for an output.

In addition, I provide two demo scripts that estimate optic flow for the same
sequence using all applicable methods. One script "DemoEstimateOpticFlowSingleLayer"
runs methods estimating a single layer motion. Another script 
"DemoEstiamteOpticFlowMultipleLayers" runs methods estimating multiple layer of
motion (here two).

Please report problems, bugs, or suggestions to
florian.raudies__at__gmail__dot__com (Replace __at__ by @ and __dot__ by .).


If you use any of the methods or code, please cite the article on scholarpedia
http://www.scholarpedia.org/article/Optic_flow.


Copyright (C) 2013 Florian Raudies, 05/16/2013, Boston University.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.