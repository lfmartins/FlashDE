package ode {
	
	/**
	* Implements the Cash-Karp method for solving systems of differential equations.
	* <p>
	* References:
	* <ul>
	* <li>J. R. Cash, A. H. Karp. "A variable order Runge-Kutta method for 
	* initial value problems with rapidly varying right-hand sides", 
	* <em>ACM Transactions on Mathematical Software</em> 16: 201-222, 1990.</li>
	* <li><a href="http://en.wikipedia.org/wiki/Cash-Karp">
	* http://en.wikipedia.org/wiki/Cash-Karp</a></li>
	* </ul>
	* </p>
	*/
	public class CashKarpSolver extends ODESolver {
		private static const  C2:Number = 1/5;
		private static const A21:Number = 1/5;
		private static const  C3:Number = 3/10;
		private static const A31:Number = 3/40;
		private static const A32:Number = 9/40;
		private static const  C4:Number = 3/5;
		private static const A41:Number = 3/10;
		private static const A42:Number = -9/10;
		private static const A43:Number = 6/5;
		private static const  C5:Number = 1; 
		private static const A51:Number = -11/54;
		private static const A52:Number = 5/2;
		private static const A53:Number = -70/27;
		private static const A54:Number = 35/27;
		private static const  C6:Number = 7/8;
		private static const A61:Number = 1631/55296;
		private static const A62:Number = 175/512;
		private static const A63:Number = 575/13824;
		private static const A64:Number = 44275/110592;
		private static const A65:Number = 253/4096;
		private static const B41:Number = 2825/27648;
		private static const B42:Number = 0;
		private static const B43:Number = 18575/48384;
		private static const B44:Number = 13525/55296;
		private static const B45:Number = 277/14336;
		private static const B46:Number = 1/4;
		private static const B51:Number = 37/378;
		private static const B52:Number = 0;
		private static const B53:Number = 250/621;
		private static const B54:Number = 125/594;
		private static const B55:Number = 0;
		private static const B56:Number = 512/1771;
		private static const  E1:Number = B51 - B41;
		private static const  E2:Number = B52 - B42;
		private static const  E3:Number = B53 - B43;
		private static const  E4:Number = B54 - B44;
		private static const  E5:Number = B55 - B45;
		private static const  E6:Number = B56 - B46;
		
		private var k2:Vector.<Number>;
		private var k3:Vector.<Number>;
		private var k4:Vector.<Number>;
		private var k5:Vector.<Number>;
		private var k6:Vector.<Number>;
		
		/**
		* Creates a solver that implements the Cash-Karp method.
		* @param system The system to be solved.
		* @param x0 Initial conditions.
		* @param t0 Initial time.
		*/
		public function CashKarpSolver(system:ODESystem=null, x0:Vector.<Number>=null, t0:Number=0) {
			super(system, x0, t0);
			derivativesPerStep = 5;
			hasErrorEstimate = true;
		}
		
		/**
		* Implements a step of the Cash-Karp method.
		*/
		protected override function step(stepsize:Number):void {			
			var t:Number = currentT + C2*stepsize;
			var tvec:Vector.<Number> = new Vector.<Number>(system.dimension);
			for (var i:int=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] + stepsize*A21*currentDX[i]
			}
			k2 = system.getDerivatives(tvec, t);

			t = currentT + C3*stepsize;
			for (i=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] + stepsize*(A31*currentDX[i] + A32*k2[i]);
			}
			k3 = system.getDerivatives(tvec, t);

			t = currentT + C4*stepsize;
			for (i=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] + stepsize*(A41*currentDX[i] + A42*k2[i] + A43*k3[i]);
			}
			k4 = system.getDerivatives(tvec, t);
			
			t = currentT + C5*stepsize;
			for (i=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] +  stepsize*(A51*currentDX[i] + A52*k2[i] + A53*k3[i] + A54*k4[i]);
			}
			k5 = system.getDerivatives(tvec, t);
			
			t = currentT + C6*stepsize;
			for (i=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] + stepsize*(A61*currentDX[i] + A62*k2[i] + A63*k3[i]  +  A64*k4[i] + A65*k5[i]);
			}
			k6 = system.getDerivatives(tvec, t);
			
			newT = currentT + stepsize;
			for (i=0; i<system.dimension; i++) {
				newX[i] = currentX[i] + stepsize*(B51*currentDX[i] +  B53*k3[i] + B54*k4[i] + B56*k6[i]);
				errorX[i] = stepsize*(E1*currentDX[i] + E3*k3[i] + E4*k4[i] + E5*k5[i] + E6*k6[i]);
			}
		}
	}
}