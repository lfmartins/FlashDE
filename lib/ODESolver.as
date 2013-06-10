package ode {
	/*
	 Class: ODESolver
	
	 "Abstract" class representing a generic
	 numerical method for systems of ODEs. To solve a system, use one
	 of the subclasses, and call either the function <solveFixed> or
	 <solveAdaptive>
	 The package has several prebuilt solvers, and users that only
	 want to solve systems of ODEs don't need to subclass this class.
	 
	 Developers that wish to define their own numerical methods must override 
	 the method <step>. 
	 
	 See also:
	 
	 - <step>
	 - <solveFixed>
	 - <solveAdaptive>
	*/
	public class ODESolver {
		private const SAFETY:Number  = 0.9;
		private const PGROW:Number   = -0.2;
		private const PSHRINK:Number = -0.25;
		private const ERRCON:Number  = Math.pow(5.0/SAFETY, 1.0/PGROW);
		private const TINY:Number    = 1.0E-30;

		private var derivativeEvaluations:uint = 0;      // Derivative evaluations in the last call of <solveFixed> or <solveAdaptive>;
		private var tolerance_ = 1E-5;                   // Bound on local error estimate
		private var stepsize_ = 0.1;                     // Time stepsize. For <solveAdaptive>, the initial step size.
		private var maxSteps_ = 10000;                   // Maximum number of steps taken in one call of <solveAdaptive>
		private var minStepsize_ = 1E-6;                 // Minimum time stepsize to be used by <solveAdaptive>
		protected var system:ODESystem;		             // The system of ODEs being solved.
		
		/*
		 Variables: 
		 
		 Protected variables important for subclassing.
		 
		 The following protected variables are important for subclasses
		 that define concrete approximation methods.
		 
		 currentT - Current value of the time variable.
		 currentX - Current values of the state variables.
		 currentDX -  Current values of the derivatives of the state variables.
		 newT - Value of the time variable after a step is performed.
		 newX - Values of the state variables after a step is performed.
		 errorX - Estimate of the local error in a step.
		 derivativesPerStep  - Number of derivative evaluations per step of the method.
		 hasErrorEstimate - Indicates if solver provides a local error estimate.
		*/
		protected var currentT:Number;                   
		protected var currentX:Vector.<Number>;          
		protected var currentDX:Vector.<Number>;        
		protected var newT:Number;                       
		protected var newX:Vector.<Number>;              
		protected var errorX:Vector.<Number>;            
		protected var derivativesPerStep:uint = 0;       
		protected var hasErrorEstimate:Boolean = false;  

		
		/**
		 Sets bound on local error estimate used by <solveAdaptive>.
		*/
		public final function set tolerance(tolerance:Number) {
			if (tolerance > 0.0) {
				tolerance_ = tolerance;
			} else {
				throw new ODEException("tolerance must be positive");
			}
		}
		
		/**
		 Gets bound on local error estimate used by <solveAdaptive>.
		*/
		public final function get tolerance():Number {
			return tolerance_;
		}
		
		/**
		 Sets time stepsize for <solveFixed>, or initial time stepsize for
		 <solveAdaptive>. The stepsize can be negative.
		*/
		public final function set stepsize(stepsize:Number) {
			if (stepsize_ != 0.0) {
				stepsize_ = stepsize;
			} else {
				throw new ODEException("stepsize cannot be zero");
			}
		}
		
		/**
		 Gets time stepsize for <solveFixed>, or initial time stepsize for
		 <solveAdaptive>
		*/
		public final function get stepsize():Number {
			return stepsize_;
		}

		/**
		 Sets maximum number of time steps taken by <solveAdaptive>.
		*/
		public final function set maxSteps(maxSteps:uint) {
			if (maxSteps > 0) {
				maxSteps_ = maxSteps;
			} else {
				throw new ODEException("maxSteps must be positive");
			}
		}
		
		/**
		 Gets maximum number of time steps taken by <solveAdaptive>.
		*/
		public final function get maxSteps():uint {
			return maxSteps_;
		}
		
		/**
		 Sets minimum stepsize used by <solveAdaptive>.
		*/
		public final function set minStepsize(minStepsize:Number) {
			if (minStepsize > 0.0) {
				minStepsize_ = minStepsize;
			} else {
				throw new ODEException("minStepSize must be positive");
			}
		}

		/**
		 Gets minimum stepsize used by <solveAdaptive>.
		*/
		public final function get minStepsize():Number {
			return minStepsize_;
		}

		/*
		 Constructor: ODESolver
		 
		 Creates an object representing a numerical method for systems of differential equations.
		 Objects of type <ODESolver> should not be constructed directly, since
		 they cannot be used to solve a system of ODEs. Use instead one of the predefined 
		 subclasses that specify a concrete numerical method.

		 This class should only be extended by programmers that wish to implement their own numerical
		 methods. Currently, only single-step methods, such as the Runge-Kutta approximations, are supported.
		 
		 Classes extending <ODESolver> must define a constructor 
		 adhering to the following conventions:
		 
		 - The constructor of the subclass must have the same signature as
		   this constuctor
		 - The constructor must call
		 
		 >super(system, x0, t0)
		
		 - The protected variable *derivativesPerStep*, if used, should be set to the
		   number of derivative evaluations performed in a single step of the method. This
		   is not necessary, unless one wishes to evaluate the performance of the method.
		
		 - The protected variable *hasErrorEstimate* should be set to 
		 *true* for numerical methods that provide a local error estimate. If
		 *hasErrorEstimate* is *false* (the default), the 
		 method <solveAdaptive> cannot be used, and throws an exception.
		 
		 Subclasses must also override the function <step>.
		 
		 Parameters:
		 
		 system - The system of differential equations to be solved.
		 x0 - The initial conditions.
		 t0 - The initial time.
		 
		 See also:
		 - <EulerSolver>
		 - <RKClassicSolver>
		 - <RKFehlbergSolver>
		 - <DormandPrinceSolver>
		 - <CashKarpSolver>
		 - <step>
		*/
		public function ODESolver(system:ODESystem=null, x0:Vector.<Number>=null, t0:Number=0) {
			if (system!=null) {
				setSystem(system,x0,t0)
			} else {
				throw new ODEException("System of differential equations cannot be null");
			}
		}

		/* Function: step
		 
		 Function defining one step of the solver.
		 
		 This function is not supposed to be directly called by users that simply want
		 to solve a given system. For this purpose, call the methods <solveFixed> or 
		 <solveAdaptive> for a subclass that implements a concrete solver.
		 
		 Subclasses defining concrete solvers must override this method, to define
		 computation of a single step of the numerical approximation. The function <step>
		 will then be called by <solveFixed> or <solveAdaptive> to compute approximate
		 solustions over an interval.
		 
		 
		 The function <step> has only an argument, the stepsize of the current step, and does not 
		 return a value.
		 For efficiency reasons, it interfaces with the solver through the
		 following protected variables:
		 
		 - currentT: the current value of the time variable.
		 - currentX: vector with the current value of the state variables.
		 - currentDX: vector with the current value of the derivatives of the state 
		 variabless.
		 - newT: the value of the time variable after the step is computed.
		 - newX: vector sith the value of the state variables after the step is 
		 computed.
		 - errorX: vector with the estimated errors in the approximation.
		 
		 
		 The variables *currentT*, *currentX* and 
		 *currentDX* are interpreted as 
		 "input variables", and should not be modified by <step>. The vector *currentDX* is guaranteed to 
		 contain the values of the derivatives of the state variables at the time *currentT*
		 and state *currentX*.
		 (This is ensured by methods <solveFixed> and <solveAdaptive>.)
		 Use of *currentDX* to access derivatives
		 is optional, but results in fewer evaluations in <solveAdaptive>;

		 The variables *newT*, *newX* and *errorX*
		 are considered "output variables". Upon exit, <step> 
		 must set *newT* and *newX*. Optionally, it 
		 may also set *errorX*, in case the numerical method has a local error estimate.
		 to be an estimate of the local error. In this case, the object's *hasErrorEstimate*
		 should be set to true in the constructor. This will allow the use of <solveAdaptive>
		 to find approximations based on an adaptive timestep.


		 Summarizing, these guidelines 
		 must be followed:
		 
		 - The variables *currentT*, *currentX* and *currentDX* should be 
		 considered "read-only".
		 - Upon return, *newT* and *newX* should be set to the value of the approximate solution
		   after a single timestep.
		 - For methods that have a local error estimate, *errorX* must also be set.
		 
		 Parameters:
		 
		 stepsize - the current time step
		 
		 Throws:
		 
		 <ODEException> if the method is called for an object of the class <ODESolver>, or
		 a subclass that does not override this function.
		*/
		protected function step(stepsize:Number):void {
			throw new ODEException("Must override function step() in subclasses");
		}
		
		/* Function: hasError
		
		 Returns true if the numerical method provides a local error estimate.

		 The function <solveAdaptive> and only be used with numerical methods that  provide
		 a local error estimate.
		
		 Returns: 
		 
		 *true* if the numerical method provides a local error estimate.
		*/
		public final function hasError():Boolean {
			return hasErrorEstimate;
		}
		
		/*
		 Function: getSystem
		 
		 Returns a reference to the <ODESystem> object associated
		 to this system.
		 
		 Do not use the returned reference change the system
		 associated to the solver. Use the method <setSystem> instead,
		 which performs proper initialization.
		
		 The main purpose of this method is to allow easy access to the parameters
		 of the system being solved. This method can be used to set and retrieve 
		 the parameters, as
		 in the examples:
		 
		 (start code)
		 	pars = solver.getSystem().parameters;
		  	solver.getSystem().paramters = {a:1.0, b:-2.0, c:0.0};
		 (end code)
		 
		 Returns:
		 
		 The system of differential equations associated to this solver object.
		 
		*/
		public final function getSystem():ODESystem {
			return system;
		}
		
		/*
		 Function: setSystem
		 
		 Sets the system of differential equations to be solved.
		 
		 Use this method to change the system of differential equations being solved. An optional
		 new initial condition can be defined too. If all that is needed is to change the
		 initial condition, this can be done by using method <setInitialConditions>
		 
		 Parameters:
		 
		  system - The system of differential equations to be solved.
		  x0 - Initial state vector.
		  t0 - Initial time.
		  
		  Throws:
		  
		  <ODEException> if system is *null*.
		  
		  See also:
		  
		  - <setInitialCondition>
		*/
		public final function setSystem(system : ODESystem, x0:Vector.<Number>=null, t0:Number=0):void {
			if (system==null) {
				throw new ODEException("System cannot be null");
			}
			this.system = system;
			currentT = t0;
			if (x0==null) {
				currentX = new Vector.<Number>(system.dimension);
				currentDX = new Vector.<Number>(system.dimension);
			} else {
				currentX = x0.slice();
				currentDX = system.getDerivatives(currentX, currentT);
			}
			newT = NaN;
			newX = new Vector.<Number>(system.dimension);
			errorX = new Vector.<Number>(system.dimension);
		}
		
		/*
		 Function: setInitialCondition
		
		 Sets the initial conditions for the system to be solved.
		 
		 Parameters:
		 
		 x0 - Initial state vector.
		 t0 - Initial time.
		*/
		public final function setInitialCondition(x0:Vector.<Number>, t0:Number):void {
			currentT = t0;
			currentX = x0.slice();
			currentDX = system.getDerivatives(currentX, currentT);
		}
		
		/*
		 Function: getCurrentT
		 
		 Returns the current value of the time variable for the system being solved.
		
		 After a solution method is called, this returns the last time at which
		 the solution was computed. This is true even if the numerical
		 approximation did not reach the final time specified in <solveAdaptive>.
		
		 Returns:
		 
		 The current value of the time variable.
		*/
		public final function getCurrentT():Number {
			return currentT;
		}

		/*
		 Function: getCurrentX
		 
		 Current value of the state variables for the system being solved.
		
		 After a solution method is called, this returns a vector with the 
		 values of of the state variables for the solution computed. 
		 This is true even if the numerical
		 approximation did not reach the final time specified in <solveAdaptive>.
		 Notice that one of the reasons that <solveAdaptive> quits is
		 loss of precision, so the values returned by this method may
		 not be reliable.
		
		 Returns: 
		 
		 The current values of the state variables.
		*/		
		public final function getCurrentX():Vector.<Number> {
			return currentX.slice();
		}

		/*
		 Function: getEvaluations
		 
		 Returns number of derivative evaluations performed in the last call of <solveFixed> or
		 <solveAdaptive>.
		 
		 This method can be used for performance measures for the numerical approximation. Ir
		 returns a count of the number of derivative evaluations in the last solution
		 computed. For the count to be accurate, the subclass defining the numerical approximation
		 must correctly set the protected varuable *derivativesPerStep*.
		 
		 Returns:
		 
		 The number of derivative evaluations in the last solution that was computed.
		*/
		public final function getEvaluations():uint {
			return derivativeEvaluations;
		}
				
		/*
		 Function: getSolutionAtPoints
		 
		 Computes the solution for a sequence of points.
		 
		 Computes an approximate solution for a sequence of points over an interval. 
		 The times are equally spaced by default. 
		 However, if maxChange os not 0, the solution attampts to add extra time values so that the
		 maximum change in the solution from one time to the next is less than maxChange. This may not be 
		 possible in some cases.
		 
		 Notice that this function does not change the values of *tolerance*, *maxSteps* and *minStepSize*
		 set by the user. These parameters should be set to appropriate values before calling this function.
		 
		 All orderings of the time values t0, t1, t2 are allowed. The computed solution will contain times
		 from t1 to t2, in this order.
		 
		 Parameters:
		 x0 - Initial values of state variables.
		 t0 - Initial time.
		 t1 - Specifies start of interval for solution.
		 t2 - Specifies end of interval for solution.
		 timeStep - The specified time step.
		 maxPoints - The maximum number of points to be computed to each side of t0. If zero, impose no limit.
		 maxChange - The maximum allowed change in the solution in a single timestep
		 minStep - Smallest step from one time value to the next.
		 adaptive - Uses the adaptive solver if true.
		 
		 Return:
		 true if computation finished successfully.
		 */
		public final function getSolutionAtPoints(x0:Vector.<Number>, t0:Number,
											   t1:Number, t2:Number, timeStep:Number,
											   xmin:Number=-1/0, xmax:Number=1/0, maxPoints:uint=0,
											   maxChange:Number=Infinity, minStep:Number=TINY, adaptive:Boolean=true):Array {
			if (maxChange <= 0) {
				throw new ODEException("maxChange must be positive");
			}

			timeStep = Math.abs(timeStep);
			if (timeStep == 0) {
				throw new ODEException("timeStep must be nonzero")
			}

			var xvalues:Array = new Array();
			var tvalues:Array = new Array();
			var npoints:uint = 0
			
			function extendSolution(tfinal:Number):Boolean {
				var lastT:Number;
				var lastX:Vector.<Number>;
				var dt:Number;
				var t:Number;
				var dist:Number;
				var xstart:Vector.<Number> = xvalues[xvalues.length - 1].slice();
				var tstart:Number = tvalues[tvalues.length - 1];
				
				var tstep:Number = timeStep;
				if (tfinal < tstart) {
					tstep = -tstep;
				}
				
				setInitialCondition(xstart,tstart);
				if (currentT == tfinal) {
					return true;
				}
				while(maxPoints == 0 || npoints <= maxPoints) {
					dt = tstep;
					lastX = currentX.slice();
					lastT = currentT;
					while (true) {
						t = currentT + dt;
						if (dt*(t-tfinal) > 0) {
							t = tfinal;
						}
						try {
							if (adaptive) {
								solveAdaptive(t);
							} else {
								solveFixed(t);
							}
						} catch (ODEException) {
							return false;
						}
						dist = supNorm(lastX, currentX)
						if (dist < maxChange) {
							break;                     
						} else {
							dt = 0.5 * dt;
							if (dt >= minStep) {
								setInitialCondition(lastX, lastT);
							} else {
								break;
							}
						}
					}
					xvalues.push(currentX.slice());
					tvalues.push(currentT);
					npoints++;
					if (dt*(currentT-tfinal) >= 0) {
						return true;
					}
					for (var i:int=0; i<currentX.length; i++) {
						if (currentX[i]<xmin || currentX[i]>xmax) {
							return false;
						}
					}
				}
				return false;
			}
			
			// We want to handle all possible orderings of t0, t1 and t2, so there are all these
			// cases. If we were lazy, we could just compute the solution at t1 and then go from
			// t1 to t2, but this is inefficient for some orderings, such as t1<t0<t2. And we are
			// not lazy.
			if ((t1 <= t0 && t0 <= t2) || (t2 <= t0 && t0 <= t1) ) {
				xvalues.push(x0.slice());
				tvalues.push(t0);
				extendSolution(t1);
				xvalues.reverse();
				tvalues.reverse();
				extendSolution(t2);
			} else if ( (t0 <= t1 && t1 <= t2) || (t2 <= t1 && t1 <= t0) ) {
				setInitialCondition(x0, t0);
				if (adaptive) {
					solveAdaptive(t1);
				}
				xvalues.push(currentX.slice());
				tvalues.push(currentT);
				extendSolution(t2);
			} else {
				setInitialCondition(x0, t0);
				if (adaptive) {
					solveAdaptive(t2);
				}
				xvalues.push(currentX.slice());
				tvalues.push(currentT);
				extendSolution(t1);
				xvalues.reverse();
				tvalues.reverse();
			}
			return [xvalues, tvalues];	   
		}
		
		private function supNorm(x1:Vector.<Number>, x2:Vector.<Number>):Number {
			var retval:Number = -Infinity;
			for (var i:uint=0; i<Math.min(x1.length, x2.length); i++) {
				 var d = Math.abs(x1[i]-x2[i]);
				 if (d > retval) {
					 retval = d
				 }
			}
			return retval;
		}
		
		/* Function: solveFixed
		 
		 Computes an approximate solution using a fixed stepsize.
		
		 Approximates the solution at a given time *t*, using a fixed stepsize.
		 Before this function is called, an initial condition must be set, either
		 in the constructor of the system being solved, or by calling
		 <setInitialCondition>. It is also necessary to set the *stepsize*.
		
		 Parameters:
		 t - The final time at which the solution is to be computed.
		 
		 Returns:
		 
		 Approximation to solution *x* of the system at time *t*, in a vector of
		 numbers.
		 
		 Throws: 
		  <ODEException> if the parameters are inconsistent an error occurs in the solution.
		*/
		public final function solveFixed(t:Number):Vector.<Number> {
			if (system==null) {
				throw new ODEException("Must define system of ODEs first");
			}
			if ((t-currentT)*stepsize <= 0) {
				stepsize = -stepsize;
			}
			
			derivativeEvaluations = 0;
			
			while ( (t - currentT) * stepsize > 0) {
				step(stepsize);
				currentT = newT;
				currentX = newX.slice();
				currentDX = system.getDerivatives(currentX, currentT);
				derivativeEvaluations += derivativesPerStep + 1;
			}
			step(t - currentT);
			currentT = newT;
			currentX = newX.slice();
			currentDX = system.getDerivatives(currentX, currentT);
			derivativeEvaluations += 1;
			
			return currentX.slice();
		}
		
		/*
		 Function: solveAdaptive
		 
		 Computes an approximate solution using a variable stepsize.
		 
		 Uses adaptive stepsize control to compute an approximate solution at
		 time *t*. Can only be used with methods that provide a local
		 error estimate, for example, Runge-Kutta methods with embedded lower order
		 approximations.
		 
		 
		 Before this function is called, an initial condition must be set, either
		 in the constructor of the system being solved, or by calling
		 <setInitialCondition>. It is also necessary to set the initial *stepsize*.
		 
		 The quality of the solution can also by setting
		 the following variables before the method is called:
		 
		 tolerance - A bound for the local error in the approximation.
		 maxSteps - The maximum number of steps to be computed.
		 minStepSize - The minimum allowed stepsize.
		 
		 Parameters:
		 
		 t - The final time at which the solution is to be computed.

		 Returns:
		 
		 Approximation to solution *x* of the system at time *t*, in a vector of
		 numbers.


		 Throws:
		 
		 <ODEException> If any of the following happens:
		 
		 - No system was defined. Before calling <solveAdaptive>, a system of differential equations
		 must be defined, either in the constructor, or by calling <setSystem>.
		 - The method does not provide a local error estimate.
		 - There is a stepsize underflow.
		 - The maximum number of iterations is exceeded.
		 - The stepsize becomes smaller than <minStepsize>
		 
		 In case an exception is thrown, information about the last approximation computed can
		 be obtained through <getCurrentT> and <getCurrentX>. However,
		 this solution may not be computed at the requested value of *t*, and
		 may not have the desired accuracy.
		 		
		References:
		
		 - W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T. Vetterling,
		 Numerical Recipes in C: The Art of Scientific Computing,
		 2ed, Cambridge Press, Cambridge, 1992.
		 
		 Available online from: 
		 <http://www.nrbook.com/nr3/>
		 
		*/
		public final function solveAdaptive(t:Number):Vector.<Number> {
			if (system==null) {
				throw new ODEException("Must define system of ODEs first");
			}
			
			if (!hasErrorEstimate) {
				throw new ODEException("Method does not provide error estimate");
			}

			if (stepsize==0) {
				stepsize = t - currentT;
			}
			
			var tdelta:Number = t - currentT;
			var stemp:Number;
			var snext:Number;
			var savedStepsize:Number;
			var restoreStepsize:Boolean = false;
			
			if (tdelta*stepsize <= 0) {
				stepsize = -stepsize
			}
			
			derivativeEvaluations = 0;
			
			var adjstep:Number;
			for (var count:int=0; count<maxSteps; count++) {
				adjstep = stepsize;
				if (tdelta*(stepsize-tdelta) > 0) {
					adjstep = tdelta;
				}
				
				var ccount:int = 0;
				while(true) {					
					step(adjstep);
					derivativeEvaluations += derivativesPerStep;
					
					var errmax:Number = 0.0;
					for (var i:int=0; i<system.dimension; i++) {
						errmax = Math.max(errmax,
										  Math.abs(errorX[i]) / 
										  (Math.abs(currentX[i]) + Math.abs(adjstep*currentDX[i]) + TINY));
					}
					errmax /= tolerance;
					if (errmax < 1.0)
						break;
					stemp = SAFETY * adjstep * Math.pow(errmax, PSHRINK);
					adjstep = (adjstep > 0) ? Math.max(stemp, 0.1*adjstep) : Math.min(stemp, 0.1*adjstep);
					
					if (currentT + adjstep == currentT) {
						throw new ODEException("Stepsize underflow");
					}
				}
				
				currentT = newT;
				currentX = newX.slice();
				currentDX = system.getDerivatives(currentX, currentT);
				derivativeEvaluations += 1;
				
				if (errmax > ERRCON) {
					stepsize = SAFETY * stepsize * Math.pow(errmax, PGROW);
				} else {
					stepsize *= 5.0
				}

				tdelta = t - currentT;
				if (Math.abs(tdelta) <= TINY) {
					break;
				}
				
				stepsize = adjstep;
				if (Math.abs(stepsize) < minStepsize) {
					throw new ODEException("Stepsize too small: " + stepsize);
				}				
			}
			
			if (count >= maxSteps) {
				throw new ODEException("Maximum number of iterations exceeded");
			}
			
			if (restoreStepsize) {
				stepsize = savedStepsize;
			}
			return currentX;
		}
	}
}



