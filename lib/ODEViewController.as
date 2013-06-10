package  ode {
	
	public class ODEViewController {
		private var system = null;
		private var solver = new CashKarpSolver();

		public function ODEViewController(system:ODESystem, solver:ODESolver=null) {
			// constructor code
			this.system = system;
			if (solver!=null) {
				this.solver = solver;
				this.solver.setSystem(system);
			}
		}

	}
	
}
