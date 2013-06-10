package ode {
	/*
	Class: ODESystem
	
	Object representing a system of ordinary differential equations.
  	
	A system of differential equations is defined by a function that computes
	the derivatives of the state variables, given a time and vector of state 
	variables. The derivatives can also, optionally, depend on given parameters. 
	*/
	public class ODESystem {
		private var dimension_:uint;
		private var derivatives:Function;
		private var parameters_:Object = new Object();
		
		/*
		 Constructor: ODESystem
		
		 Creates an object representing a system of ordinary differentials equations.
		 
		 A system of ODEs is defined by the function that computes the derivatives, 
		 the dimension of the system, and (optionally) parameters.
		 
		 The function that computes the derivatives must have 
		 the following signature:
		 > function fsys(x:Vector.<Number>, t:Number):Vector.<Number>
		 
		 The following rules 
		 should be observed:
		 
		 - The first argument of *fsys* must be a vector with length equal to
		 the dimension of the system.
		 - The return value of *fsys* must be a vector with length equal to the
		 dimentsion of the system.
		 - The *parameters* must be an object whose fields are numbers, for example:
		 
		 > {alpha:0.5, beta:-1.2}
		 
		 In this example, the function *fsys* would use *parameters.alpha*
		 and *parameters.beta* to refer to the parameters. It is not possible 
		 to use "deep" objects as parameters, such as:
		 
		 > {v:[2.3,4.5]}     // ERROR.
		
		 Parameters:
		
		 fsys - Function defining the system.
		
		 dimension - Dimension of the system.
		
		 parameters - Parameters of the system.
		 
		 Example:
		 
		 The code below defines a system representing 
		 a forced harmonic oscillator:
		 
		 \[\frac{dx_0}{dt} = x_1, \quad \frac{dx_1}{dt}=-kx_0-cx_1+A \sin(\omega t).\]
		 
		 Notice how parameters are passed in an *Object*. Then, the code calls the function
		 *fsys* to compute the value of the field at a given point and time. This example does
		 not solve the system.
		 
		(code)
		// File: ForcedOscillator.as
		package  {
			import flash.display.MovieClip;
			import ode.ODESystem;
	
			public class ForcedOscillator01 extends MovieClip {
				private var system:ODESystem;

				public function ForcedOscillator01()  {
					system = new ODESystem(fsys, 2, {k:1.0, c:0.0, A:2.0, omega:Math.PI});
			
					var xvec:Vector.<Number> = Vector.<Number>([1.0,2.0]);
					var t:Number = 1.0;
					var dxvec:Vector.<Number> = system.getDerivatives(xvec, t);
			
					trace( "The field at the point [" + xvec +
						   "] at time t=" + t + " is [" + dxvec + "]" );
				}

				private function fsys(x:Vector.<Number>, t:Number) {
					var p:Object = system.parameters
					return Vector.<Number>( [x[1], 
											 - p.k*x[0] - p.c*x[1] + p.A * Math.sin(p.omega*t)]);
			}
			}
		}
		(code)
		*/
		public function ODESystem(fsys:Function, dimension:int, parameters:Object=null) {
			this.dimension_ = dimension;
			this.derivatives = fsys;
			this.parameters = parameters;
		}
		
		/*
		  Function: getDerivatives
		  
		  Computes the derivatives of the state vector at a given time.
		 
		  Returns a vector with the derivatives of the state variable at a given time
		  (in other words, it computes the "right hand side" of the system differential
		  equations.
		 
		  Parameters:
		
		  x - State vector at which derivatives are computed.
		  t - Time at which derivatives are computed.
		
		  Returns:
		 
		  Vector with the derivatives of the state vector *x*, 
		  at time *t*.
		
		  Throws:
		 
		  An <ODEException> if an error occurs in computing the derivatives 
		*/
		public function getDerivatives(x:Vector.<Number>, t:Number):Vector.<Number> {
			if (x.length != dimension) {
				throw new ODEException("Vector x has wrong dimension");
			}
			try {
				var dx:Vector.<Number> = derivatives(x, t);
			} catch (e:ArgumentError) {
				throw new ODEException("Function derivatives() must have two arguments");
			} catch (e:RangeError) {
				throw new ODEException("Reference to x[i] with i out of range");
			} catch (e:TypeError) {
				throw new ODEException("Function derivatives() has argument or " + 
									   "return value with wrong type");
			} catch (e:ReferenceError) {
				throw new ODEException("Reference to undefined property. " +
									   "Check parameters in function derivatives()");
			}
			
			if (dx.length != dimension) {
				throw new ODEException("Return value of function derivatives() has " +  
									   "wrong dimension");
			}
			
			for (var i:int=0; i<dimension; i++) {
				if ( !isFinite(dx[i]) ) {
					throw new ODEException("Derivatives are infinite or NaN");
				}
			}
			
			return dx;
		}
		
		/*
		 Function: set parameters
		
		 Sets the parameters for this system.
		 
		 Sets the parameters that will, subsequently, be used for computations with the system.
		 The *parameters* must have the same fields used in the function *fsys*
		 used to define the system.
		
		 This function does a "shallow copy" of the fields in the *parameters* object. Since
		 the fields of the *parameters* object should be of type Number, this is never a
		 problem.
		
		 Parameters:
		
		 parameters - The new value of the parameters.
		*/
		public function set parameters(parameters:Object):void {
			for (var v:* in parameters) {
				if (! (parameters[v] is Number) ) {
					throw new ODEException( "Parameters must be of type number" );
				}
											   
			}
			parameters_ = new Object();
			// Do a shallow copy for parameters
			for (v in parameters) {
				parameters_[v] = parameters[v]
			}
		}
		
		/*
		 Function: get parameters
		
		 Returns a shallow copy of the parameters being used by this system.
		
		 Returns:
		 
		 The parameters currently in use by the system of differential equations.
		*/
		public function get parameters():Object {
			var retval:Object = new Object();
			
			for (var v:* in parameters_) {
				retval[v] = parameters_[v]
			}
			return retval;
		}
		
		/*
		 Function: get dimension
		 
		 Returns the dimension of the system.
		 
		 It is not possible to reset the
		 dimension of a system once it has been defined.
		
		 Returns:
		 
		 The dimension of the system. 
		*/
		public function get dimension():uint {
			return dimension_;
		}
	}
}

