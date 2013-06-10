package ode {
	import flashandmath.as3.parsers.*
	
	/*
	 Class: SymbolicODESystem
	
	 Object representing a system of ordinary differential equations
	 defined in terms of symbolic expressions. 
	 This class requires the Flash and Math library, avaliable at:
	 <http://www.flashandmath.com/>. Mathematical formulas use the
	 syntax defined by the Flash&Math parser. See the references below.

	References:
	
	 - D. Ensley, B.Kaskosz. Flash and Math Applets, Learn by Example: Programming 
	   in ActionScript 3 for Mathematics and Science Teaching and Learning (2009).</li>
	 - *flashandmath* package documentation: 
	 <http://www.flashandmath.com/basic/simplegraph/flashandmath_guide.pdf>
	*/
	public class SymbolicODESystem extends ODESystem {
		private var parser:MathParser;
		private var compiled:Array;
		private var dexps:Array;
		private var vnames:Array;
		private var pnames:Array;
		private var allnames:Array;
		
		/*
		 Constructor: SymbolicODESystem
		 
		 Defines a system of ordinary differential equations defined in terms of symbolic expressions.
		
		 Parameters:
		 
		 vnames - Array of strings representing the variables used in the system definition.  
		          The first string (index 0) is the independent variable, the others
				  are the state variables.
		 dexps - Array of strings representing the expressions defining the 
		         derivatives of the state variables. Include only the expression for
				 the right-hand side of the equations, following the order specified in the
				 array *vnames*. The number of expressions must be equal to the number of
				 state variables.
		 pnames - Array of strings representing the parameters used in the system definition.
		          
		 parameters - Vector of numbers used to initialize parameters. Must have the same
		              length as the array *pnames*, and the parameter values must be given
					  in the order specified by this array.
					  
		 Throws:
		  
		 <ODEException> If there is a syntax error in the symbolic expressions.
		 
		 Example:
		 
		 The code below defines a system representing 
		 a forced harmonic oscillator:
		 
		 \[\frac{dx_0}{dt} = x_1, \quad \frac{dx_1}{dt}=-kx_0-cx_1+A \sin(\omega t).\]
		(code)
		package  {
			import flash.display.MovieClip;
			import ode.SymbolicODESystem;
	
			public class ForcedOscillatorSymbolic extends MovieClip {
				var system:SymbolicODESystem;

				public function ForcedOscillatorSymbolic() {
					system = new SymbolicODESystem(['t','x','v'], 
												   ['v', '-k*x-c*v+A*sin(omega*t)'],
												   ['k', 'c', 'A', 'omega'],
												   Vector.<Number>([1.0, 0.0, 2.0, Math.PI]));
					var xvec:Vector.<Number> = Vector.<Number>([1.0,2.0]);
					var t:Number = 1.0;
					var dxvec:Vector.<Number> = system.getDerivatives(xvec, t);
					trace( "The field at the point [" + xvec +
						   "] at time t=" + t + " is [" + dxvec + "]" );
				   
				    system.parameters = {k:2.0, c:1.0, A:-3.0, omega:3*Math.PI};
					xvec = Vector.<Number>([-3.0, 5.0]);
					t = 1.25;
					dxvec = system.getDerivatives(xvec, t);
					trace( "The field at the point [" + xvec +
						   "] at time t=" + t + " is [" + dxvec + "]" );			
				}
			}
		}
		(code)
		*/
		public function SymbolicODESystem(vnames:Array, dexps:Array, pnames:Array, parameters:Vector.<Number>) {
			super(_derivatives, dexps.length);
			if (vnames.length-1 != dexps.length) {
				throw new ODEException("Array dexps must have one expression for each state variable");
			}
			if (pnames.length != parameters.length) {
				throw new ODEException("Arrays pnames and parameters must have the same length");
			}
			this.vnames = vnames.slice();
			this.pnames = pnames.slice();
			this.dexps = dexps.slice();
			allnames = vnames.concat(pnames);
			parser = new MathParser(allnames);
			compiled = new Array(dexps.length);
			for (var i:int=0; i<dexps.length; i++) {
				compiled[i] = parser.doCompile(dexps[i]);
				if (compiled[i].errorStatus==1) {
					throw new ODEException("Error in symbolic system: " + compiled[i].errorMes);
				}
			}
			var pobj = new Object();
			for (i=0; i<pnames.length; i++) {
				pobj[pnames[i]] = parameters[i];
			}
			this.parameters = pobj;
		}
		
		private function _derivatives(x:Vector.<Number>, t:Number):Vector.<Number> {
			var dx:Vector.<Number> = new Vector.<Number>(dimension);
			var values:Array = [t];
			for (var i:int=0; i<dimension; i++) {
				values.push(x[i]);
			}
			for (i=0; i<pnames.length; i++) {
				values.push(parameters[pnames[i]]);
			}
			for (i=0; i<dimension; i++) {
				dx[i] = parser.doEval(compiled[i].PolishArray,values);
			}
			return dx;
		}
	}
}
