package ode {
	
	/**
	* Implements the Dormand-Prince method for solving systems of differential equations.
	* <p>
	* References:
	* <ul>
	* <li>Dormand, J. R.; Prince, P. J., "A family of embedded Runge-Kutta formulae", 
	* <em>Journal of Computational and Applied Mathematics 6 (1): 19–26, 1980</em></li>
	* <li><a href="http://en.wikipedia.org/wiki/Dormand-Prince">
	* http://en.wikipedia.org/wiki/Dormand-Prince</a></li>
	* </ul>
	* </p>
	*/
	public class DormandPrinceSolver extends ODESolver {
		private static const  C2:Number = 1/5;
		private static const A21:Number = 1/5;
		private static const  C3:Number = 3/10;
		private static const A31:Number = 3/40;
		private static const A32:Number = 9/40;
		private static const  C4:Number = 4/5;
		private static const A41:Number = 44/45;
		private static const A42:Number = -56/15;
		private static const A43:Number = 32/9;
		private static const  C5:Number = 8/9;
		private static const A51:Number = 19372/6561;
		private static const A52:Number = -25360/2187;
		private static const A53:Number = 64448/6561;
		private static const A54:Number = -212/729;
		private static const  C6:Number = 1;
		private static const A61:Number = 9017/3168;
		private static const A62:Number = -355/33;
		private static const A63:Number = 46732/5247;
		private static const A64:Number = 49/176;;
		private static const A65:Number = -5103/18656;
		private static const  C7:Number = 1;
		private static const A71:Number = 35/384;
		private static const A72:Number = 0;
		private static const A73:Number = 500/1113;
		private static const A74:Number = 125/192;
		private static const A75:Number = -2187/6784;
		private static const A76:Number = 11/84;
		private static const B41:Number = 1951/21600;
		private static const B42:Number = 0;
		private static const B43:Number = 22642/50085;
		private static const B44:Number = 451/720;
		private static const B45:Number = -12231/42400;
		private static const B46:Number = 649/6300;
		private static const B47:Number = 1/60;
		private static const B51:Number = 35/384;
		private static const B52:Number = 0;
		private static const B53:Number = 500/1113;
		private static const B54:Number = 125/192;
		private static const B55:Number = -2187/6784;
		private static const B56:Number = 11/84;
		private static const B57:Number = 0
		private static const  E1:Number = B51 - B41;
		private static const  E2:Number = B52 - B42;
		private static const  E3:Number = B53 - B43;
		private static const  E4:Number = B54 - B44;
		private static const  E5:Number = B55 - B45;
		private static const  E6:Number = B56 - B46;
		private static const  E7:Number = B57 - B47;
		
		private var k2:Vector.<Number>;
		private var k3:Vector.<Number>;
		private var k4:Vector.<Number>;
		private var k5:Vector.<Number>;
		private var k6:Vector.<Number>;
		private var k7:Vector.<Number>;
		
		/**
		* Creates a solver that implements the  Dormand-Prince method.
		* @param system The system to be solved.
		* @param x0 Initial conditions.
		* @param t0 Initial time.
		*/
		public function DormandPrinceSolver(system:ODESystem=null, x0:Vector.<Number>=null, t0:Number=0) {
			super(system, x0, t0);
			derivativesPerStep = 6;
			hasErrorEstimate = true;
		}
		
		/**
		* Implements a step of the Dormand-Prince method.
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

			t = currentT + C7*stepsize;
			for (i=0; i<system.dimension; i++) {
				tvec[i] = currentX[i] + stepsize*(A71*currentDX[i] + A72*k2[i] + A73*k3[i]  +  A74*k4[i] + A75*k5[i] + A76*k6[i]);
			}
			k7 = system.getDerivatives(tvec, t);
			
			newT = currentT + stepsize;
			for (i=0; i<system.dimension; i++) {
				newX[i] = currentX[i] + stepsize*(B51*currentDX[i] +  B52*k2[i] +B53*k3[i] + B54*k4[i] + B55*k5[i] + B56*k6[i] + B57*k7[i]);
				errorX[i] = stepsize*(E1*currentDX[i] + E2*k2[i] +E3*k3[i] + E4*k4[i] + E5*k5[i] + E6*k6[i] + E7*k7[i]);
			}
		}
	}
}