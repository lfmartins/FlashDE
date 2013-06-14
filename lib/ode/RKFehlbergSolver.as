package ode {
	
	/**
	* Implements the Runge-Kutta-Fehlberg method for solving systems of differential equations.
	* <p>
	* References:
	* <ul>
	* <li>
	* W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T. Vetterling,
	* <em>Numerical Recipes in C: The Art of Scientific Computing</em>,
	* 2ed, Cambridge Press, Cambridge, 1992.</li>
	* <li><a href="http://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method">
	* http://en.wikipedia.org/wiki/Runge-Kutta-Fehlberg_method</a></li>
	* </ul>
	* </p>
	*/
	public class RKFehlbergSolver extends ODESolver {
		private static const  C2:Number = 1/4;
		private static const A21:Number = 1/4;
		private static const  C3:Number = 3/8;
		private static const A31:Number = 3/32;
		private static const A32:Number = 9/32;
		private static const  C4:Number = 12/13;
		private static const A41:Number = 1932/2197;
		private static const A42:Number = -7200/2197;
		private static const A43:Number = 7296/2197;
		private static const  C5:Number = 1;
		private static const A51:Number = 439/216;
		private static const A52:Number = -8;
		private static const A53:Number = 3680/513;
		private static const A54:Number = -845/4104;
		private static const  C6:Number = 1/2;
		private static const A61:Number = -8/27;
		private static const A62:Number = 2;
		private static const A63:Number = -3544/2565;
		private static const A64:Number = 1859/4104;
		private static const A65:Number = -11/40;
		private static const B41:Number = 25/216;
		private static const B42:Number = 0;
		private static const B43:Number = 1408/2565;
		private static const B44:Number = 2197/4104;
		private static const B45:Number = -1/5;
		private static const B46:Number = 0;
		private static const B51:Number = 16/135;
		private static const B52:Number = 0;
		private static const B53:Number = 6656/12825;
		private static const B54:Number = 28561/56430;
		private static const B55:Number = -9/50;
		private static const B56:Number = 2/55;
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
		
		public function RKFehlbergSolver(system:ODESystem=null, x0:Vector.<Number>=null, t0:Number=0) {
			super(system, x0, t0);
			derivativesPerStep = 5;
			hasErrorEstimate = true;
		}
		
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
				newX[i] = currentX[i] + stepsize*(B51*currentDX[i] + B53*k3[i] + B54*k4[i] + B55*k5[i] + B56*k6[i]);
				errorX[i] = stepsize*(E1*currentDX[i] + E3*k3[i] + E4*k4[i] + E5*k5[i] + E6*k6[i]);
			}
		}
	}
}