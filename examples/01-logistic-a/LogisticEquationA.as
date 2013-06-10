package  {
	import flash.display.MovieClip;
	import ode.ODESystem;
	
	/* Class: LogisticEquationA
	
	Class for the logistic-a.fla example file.
	
	To use this example, open logistic-a.fla in Flash.
	
	This examples defines the differential equation corresponding
	to the logistic population growth model and computes the derivative
	at a given point. Solutions are not actually computed. Here is the code
	for the class:
	
	(start code)
	public class LogisticEquationA extends MovieClip {
		private var system:ODESystem;
		
		public function LogisticEquationA()  {
			system = new ODESystem(logistic_equation, 1);
			
			var xvec:Vector.<Number> = Vector.<Number>([2.0]);
			var t:Number = 1.0;
			var dxvec:Vector.<Number> = system.getDerivatives(xvec, t);
			
			trace( "The derivative dx/dt at the point x = " + xvec[0] + 
				   " and time " + t + " is " + dxvec[0]);
		}

		private function logistic_equation(x:Vector.<Number>, t:Number) {
			return new <Number>[x[0]*(1-x[0])];
		}
	}	
	(end code)
	*/
	public class LogisticEquationA extends MovieClip {
		private var system:ODESystem;
		
		public function LogisticEquationA()  {
			system = new ODESystem(logistic_equation, 1);
			
			var xvec:Vector.<Number> = Vector.<Number>([2.0]);
			var t:Number = 1.0;
			var dxvec:Vector.<Number> = system.getDerivatives(xvec, t);
			
			trace( "The derivative dx/dt at the point x = " + xvec[0] + 
				   " and time " + t + " is " + dxvec[0]);
		}

		private function logistic_equation(x:Vector.<Number>, t:Number) {
			return new <Number>[x[0]*(1-x[0])];
		}
	}
}
