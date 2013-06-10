package  {
	import flash.display.MovieClip;
	import ode.ODESystem;
	import ode.RKFehlbergSolver;
	
	/* Class: LogisticEquationB
	
	Class for the logistic-b.fla example file.
	
	To use this example, open logistic-b.fla in Flash.
	
	This examples defines the differential equation corresponding
	to the logistic population growth model and computes an approximate
	solution using the Runge-Kutta-Fehlberg solver.
	*/
	public class LogisticEquationB extends MovieClip {
		private var system:ODESystem;
		private var solver:RKFehlbergSolver;
		
		public function LogisticEquationB()  {
			system = new ODESystem(logistic_equation, 1);
			solver = new RKFehlbergSolver(system);
			
			var x0:Vector.<Number> = Vector.<Number>([2.0]);
			var t0:Number = 0.0;
			var tmax:Number = 10.0;
			
			solver.tolerance = 1E-3;
			solver.stepsize = 0.1;
			solver.maxSteps = 500;
			var solution:Array = solver.getSolutionAtPoints(x0, t0, t0, tmax, 1.0);
			trace(solution.length);
			var xvalues:Array = solution[0];
			var tvalues:Array = solution[1];
			
			for (var i:int = 0; i<tvalues.length; i++) {
				var t:Number = tvalues[i];
				var xval:Number = xvalues[i][0];
				var xexact:Number = 1.0 / (1 + (1/x0[0]-1)*Math.exp(-t));
				trace( "t = " + t + " x(t) = " + xval + " Exact = " + xexact + " Error = " + Math.abs(xval-xexact));
			}
		}

		private function logistic_equation(x:Vector.<Number>, t:Number) {
			return new <Number>[x[0]*(1-x[0])];
		}
	}
}
