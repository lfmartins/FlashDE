package ode {
	
	/*
	 Class: ODEException
	 
	 Exception class for ODEs library.
	*/
	public class ODEException extends Error {
		
		public function ODEException(message:String) {
			super(message);
			name = "ODEException";
		}
	}
}