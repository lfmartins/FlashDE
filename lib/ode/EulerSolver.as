package ode {
	
	/**
	* Implements the Euler's method for solving systems of differential equations.
	*/
	public class EulerSolver extends ODESolver {
		
		/**
		* Creates a solver that implements the  Euler's method.
		* @param system The system to be solved.
		* @param x0 Initial conditions.
		* @param t0 Initial time.
		*/
		public function EulerSolver(system:ODESystem=null, x0:Vector.<Number>=null, t0:Number=0) {
			super(system, x0, t0);
			derivativesPerStep = 0;
			hasErrorEstimate = false;
		}
		
		/**
		* Implements a step of Euler's method.
		*/
		protected override function step(stepsize:Number):void {
			newT = currentT + stepsize;
			for (var i:int=0; i<system.dimension; i++) {
				newX[i] = currentX[i] + stepsize*currentDX[i]
			}
		}
	}
}