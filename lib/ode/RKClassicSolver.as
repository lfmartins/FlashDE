package ode {
	/**
	Class: RKClassicSolver
	
	Implements the classical 4th order Runge-Kutta method for solving
	Differential equations.
	
	References:		
	 - W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T. Vetterling, Numerical Recipes in C: The Art of Scientific Computing,
	 2ed, Cambridge Press, Cambridge, 1992.
	
	 - <a href="http://en.wikipedia.org/wiki/Runge-Kutta">
	   http://en.wikipedia.org/wiki/Runge-Kutta</a></li>
	*/
	public class RKClassicSolver extends ODESolver {
		private var k2:Vector.<Number>;
		private var k3:Vector.<Number>;
		private var k4:Vector.<Number>;
		
		/**
		* Creates a solver that implements the classic 4th order Runge-Kutta method.
		* @param system The system to be solved.
		* @param x0 Initial conditions.
		* @param t0 Initial time.
		*/
		public function RKClassicSolver(system:ODESystem=null, x0:Vector.<Number>=null, t0:Number=0) {
			super(system, x0, t0);
			derivativesPerStep = 3;
			hasErrorEstimate = false;
		}

		/**
		* Implements a step of the classic 4th order Runge-Kutta method.
		*/
		protected override function step(stepsize:Number):void {
			var t:Number = currentT + stepsize/2;			
			var tvec:Vector.<Number> = new Vector.<Number>(system.dimension);
			for (var i:int=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] + stepsize * currentDX[i]/2
			}
			k2 = system.getDerivatives(tvec, t);
			
			for (i=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] + stepsize * k2[i]/2;
			}
			k3 = system.getDerivatives(tvec, t);
			
			t = currentT + stepsize;
			for (i=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] + stepsize * k3[i];
			}
			k4 = system.getDerivatives(tvec, t);
			
			newT = t;
			for (i=0; i<system.dimension; i++) {
				newX[i] = currentX[i] + stepsize * ((currentDX[i]+k4[i])/6 + (k2[i]+k3[i])/3);
			}
		}
	}
}