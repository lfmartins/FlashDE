﻿package ode {	import flash.display.GraphicsPathCommand;	import flash.display.Shape;	import flash.display.Sprite;	import flash.display.GraphicsPathCommand;		import flash.text.TextField;	import flash.text.TextFormat;	import flash.text.TextFieldType;	import flash.text.TextFieldAutoSize;	import flash.text.AntiAliasType;		import flash.utils.Timer;	import flash.events.TimerEvent;	import flash.events.MouseEvent;		import flashandmath.as3.boards.GraphingBoard;	import ode.ODESystem;	import ode.ODEException;		/*	 Class: ODEBoard	 	 An extension of GraphicsBoard specialized to plotting solutions of ODEs.	 	 This class extends the GraphicsBoard class from the Flash&Math library, adding	 functionality to plot solutions of ODEs. Currently, two kinds of plots are supported:	 time-dependent plots, and plots of two-dimensional projections of the phase space.	 	 The class has facilities for plotting solution curves and direction fields. Solutions	 are plotted and stored, so that numerical values can be accessed after computation.	*/	public class ODEBoard extends GraphingBoard {		private const TIMEDEPENDENT:uint = 0;		private const PHASESPACE:uint = 1;		private const DEFAULT_THICKNESS:uint = 2;		private const DEFAULT_COLOR:uint = 0x000000;		private const DEFAULT_ALPHA:uint = 1.0;		private const MAXCOMPONENTS:uint = 30;				private var xmin:Number;		private var xmax:Number;		private var ymin:Number;		private var ymax:Number;		private var tmin:Number;		private var tmax:Number;				private var arrows:Array;		private var arrowsX:Array;		private var arrowsY:Array;				private var dxArrows:Number;		private var dyArrows:Number;				private var system:ODESystem;		private var solver:ODESolver;				private var arrowColor:uint = 0x000000;		private var arrowShaftWidth:Number = 8;		private var arrowPointWidth:Number = 6;		private var arrowHeight:Number = 5;		private var arrowThickness:Number = 2;		private var arrowReferencePoint:int = Arrow.REF_BASE;				private var spSolutions:Sprite = new Sprite();		private var styles:Array  = new Array(MAXCOMPONENTS);				private var particles:Array = new Array();		private var particleSize:uint = 5;		private var particleColor:uint = 0x99CCFF;		private var showParticles_:Boolean = true;				private var keepSolutionData:Boolean;		private var solutionData:Array = new Array();		private var drawSolutionOnClick:Boolean = true;		private var initialConditions:Array = new Array();						private var tfMessage:TextField = new TextField();				private var plotComponents:Array;		private var plotType:uint;		private var hasArrows:Boolean;						/*		 Constructor: ODEBoard		 		 Creates an object of this class.		 		 Creates and <ODEBoard>. If a system is given, it initializes the board with a 		 the corresponding direction field. Notice that only the following kinds of plots		 support a direction field:		 - Time-dependent plots of one-dimensional systems.		 - Phase space plots of two-dimensional systems.		 		 The components to be plotted are specified in the array *plotComponents*. This array		 can be one of the following:		 		 - A time-dependent plot. This is indicated by an array of integers, each integer indexing		   one component to be plotted.		 - A two-dimensional projection of the phase space. This is indicated by an array of pairs		   of integers.		   		 If *plotComponents* is *null*, 		 the following defaults are used:		 		 - For one-dimensional systems, a time-dependent graph of the solution.		 - For two-dimensional systems, a phase space plot.		 - For systems with more than two dimesions, a time-dependent plot of all components.		 		 Parameters:		 		 - width: width in pixels.		 - height: height in pixels.		 - system: system of ODEs (cannot be *null*).		 - plotComponents: specifies components to be plotted.		 		 Throws:		 		 <ODEException> if *system* is null.		*/		public function ODEBoard(width:Number, height:Number, system:ODESystem,								 plotComponents:Array=null, keepSolutionData:Boolean=false) {			super(width, height);			if (system==null) {				throw new ODEException("system cannot be null")			}			this.system = system;			this.plotType = plotType;			spSolutions.x = 0;			spSolutions.y = 0;			addChild(spSolutions);			var smask:Shape = new Shape();			addChild(smask);			with (smask.graphics) {				beginFill(0x000000);				drawRect(0,0,width,height);			}			spSolutions.mask = smask;						var format:TextFormat = new TextFormat();			format.font="Arial";    		format.size=14;									if (plotComponents==null || plotComponents.length==0) {				if (system.dimension==1) {					plotComponents = [0];				} else if (system.dimension==2) {					plotComponents = [[0,1]];				} else {					plotComponents = new Array(system.dimension);					for (i=0; i<system.dimension; i++) {						plotComponents[i] = i;					}				}			}			if (plotComponents.length > MAXCOMPONENTS) {				throw new ODEException("Maximum number of components allowed is " + MAXCOMPONENTS); 			}						this.plotComponents = plotComponents;			for (var i:uint=0; i < MAXCOMPONENTS; i++) {				styles[i] = [DEFAULT_THICKNESS, DEFAULT_COLOR, DEFAULT_ALPHA]			}						hasArrows = false;			if (plotComponents[0] is uint) {				for (i=0; i<plotComponents.length; i++) {					if (! (plotComponents[i] is uint) || plotComponents[i]>=system.dimension){						throw new ODEException('plotComponents are invalid');					}				}				plotType = TIMEDEPENDENT;				if (system.dimension==1) {					hasArrows = true;				}			} if (plotComponents[0] is Array) {				for (i=0; i<plotComponents.length; i++) {					var plc = plotComponents[i];					if (!(plc is Array) || 						plc.length!=2 ||						!(plc[0] is uint) ||						!(plc[1] is uint) ||						plc[0] >= system.dimension ||						plc[1] >= system.dimension ||						plc[0] == plc[1]) {						throw new ODEException('plotComponents are invalid');					}				}				plotType = PHASESPACE;				if (system.dimension==2) {					hasArrows = true;				}			}						setVarsRanges(-10,10,-10,10);			setTimeRange(-10,10);			if (plotType==PHASESPACE) {				setArrowProperties(10,6,5,3,0x000000,Arrow.REF_BASE);			} else {				setArrowProperties(10,0,0,3,0x000000,Arrow.REF_MIDDLE);			}			setArrows(1, 1, -9, 9, -9, 9);			setArrowsColor(0x000000);			updateArrows(0);						solver = new CashKarpSolver(system);						addEventListener (MouseEvent.CLICK, onMouseClick);		}				public function getSystem():ODESystem {			return system;		}				public function setSystem(system:ODESystem):void {			this.system = system;			if (system==null) {				throw new ODEException("system cannot be null")			}		}				public override function setVarsRanges(xmin:Number,xmax:Number,ymin:Number,ymax:Number):void {			super.setVarsRanges(xmin,xmax,ymin,ymax);			this.xmin = xmin;			this.xmax = xmax;			this.ymin = ymin;			this.ymax = ymax;		}				public function setTimeRange(tmin:Number, tmax:Number) {			this.tmin = tmin;			this.tmax = tmax;		}				public function setArrowProperties(shaftWidth:Number, pointWidth:Number, height:Number,										   thickness:Number, color:Number, reference:uint):void {			arrowShaftWidth = shaftWidth;			arrowPointWidth = pointWidth;			arrowHeight = height;			arrowThickness = thickness;			arrowColor = color;			arrowReferencePoint = reference;		}				public function setArrows(dx:Number, dy:Number, xmin:Number=NaN, xmax:Number=NaN,								  ymin:Number=NaN, ymax:Number=NaN):void {			var i:int;			var j:int;			var xx:Number;			var yy:Number;			if (dx<=0 || dy<=0) {				throw new ODEException("dx and dy must be positive.")			}			if (!hasArrows) {				return;			}			if (arrows!=null) {				for (i=0; i<arrows.length; i++) {					for (j=0; j<arrows[i].length; j++) {						if (arrows[i][j] != null) {							removeChild(arrows[i][j]);						}					}				}			}			if (plotType == TIMEDEPENDENT) {				if (isNaN(xmin)) {xmin = this.xmin;}				if (isNaN(xmax)) {xmax = this.xmax;}				if (isNaN(ymin)) {ymin = this.ymin;}				if (isNaN(ymax)) {ymax = this.ymax;}			} else if (plotType == PHASESPACE) {				if (isNaN(xmin)) {xmin = this.tmin;}				if (isNaN(tmax)) {xmax = this.tmax;}				if (isNaN(ymin)) {ymin = this.xmin;}				if (isNaN(ymax)) {ymax = this.xmax;}			} else {				throw new ODEException("invalid plot type.");			}			if (xmin >= xmax || ymin >= ymax) {				throw new ODEException("must have xmin < xmax and ymin < ymax");			}			arrowsX = new Array();			arrowsY = new Array();			for (i=0; (xx = xmin + i*dx) <= xmax; i++) {				arrowsX.push(xx); 			}			for (j=0; (yy = ymin + j*dy) <= ymax; j++) {				arrowsY.push(yy); 			}			arrows = new Array(arrowsX.length);			for (i=0; i<arrowsX.length; i++) {				arrows[i] = new Array(arrowsY.length);				for (j=0; j<arrowsY.length; j++) {					var arrow:Arrow = new Arrow(arrowShaftWidth, arrowPointWidth, arrowHeight, 												arrowThickness, arrowColor, arrowReferencePoint);										arrow.x = xtoPix(arrowsX[i]);					arrow.y = ytoPix(arrowsY[j]);					addChildAt(arrow,getChildIndex(spSolutions)-1);										arrows[i][j] = arrow;				}			}		}				public function updateArrows(t:Number):void {			if (!hasArrows) {				return;			}			for (var i:int=0; i<arrows.length; i++) {				for (var j:int=0; j<arrows[i].length; j++) {					try {						var dxdy:Vector.<Number>;						if (plotType == PHASESPACE) {							dxdy = system.getDerivatives(Vector.<Number>([arrowsX[i],arrowsY[j]]), t);						} else {							dxdy = system.getDerivatives(Vector.<Number>([arrowsY[j]]), arrowsX[i]);						}					} catch (e:ODEException) {						continue;					}					var dx:Number;					var dy:Number;					if (plotType == PHASESPACE) {						var c:Array = plotComponents[0];						var ix:int = c[0];						var iy:int = c[1];						dx = dxdy[ix];						dy = dxdy[iy];					} else {						ix = plotComponents[0];						dx = 1.00;						dy = dxdy[ix];					}					var angle:Number;					if (dx==0) {						angle = (dy>=0) ? Math.PI/2 : -Math.PI/2;					} else { 						angle = Math.atan(dy/dx);						if (dx<0) {							angle = Math.PI+angle;						}					}					arrows[i][j].isArrow = Math.abs(dx)+Math.abs(dy)>1E-5;					arrows[i][j].rotation = -angle*180/Math.PI;				}			}		}				public function set showParticles(b:Boolean):void {			showParticles_ = b;			for (var i:int=0; i<particles.length; i++) {				particles[i].visible = b;			}		}				public function updateParticles(t:Number):void {			if (!showParticles_) {				return;			}			for (var i:int=0; i<solutionData.length; i++) {				var x:Vector.<Number> = getSolutionAt(i,t);				if (x==null) {					particles[i].visible = false;				} else {					particles[i].visible = true;					particles[i].x = xtoPix(x[0]);					particles[i].y = ytoPix(x[1]);				}			}		}				public function setArrowsVisible(b:Boolean=true):void {			if (!hasArrows) {				return;			}			for (var i:int=0; i<arrows.length; i++) {				for (var j:int=0; j<arrows[i].length; j++) {					arrows[i][j].visible = b;				}			}		}				public function setArrowsColor(color:uint):void {			if (!hasArrows) {				return;			}			arrowColor = color;			for (var i:int=0; i<arrows.length; i++) {				for (var j:int=0; j<arrows[i].length; j++) {					arrows[i][j].color = color;				}			}		}				public function setParticlesColor(color:uint):void {					}						/*			Function: setStyles						Sets the styles for solution curves.						The input for this function is an array of 3-compnent arrays, specifying the line style			for each plotted component.			The line style for each plotted component is specifyed in an array with three components, 			indicating the line's thickness, color and alpha.		*/		public function setStyles(styles:Array):void {			for (var k:uint=0; k < Math.min(MAXCOMPONENTS, styles.length); k++) {				var style:Array = styles[k];				if ( !(style is Array) || style.length != 3 || 					 !(style[0] is Number) || !(style[1] is uint) || !(style[2] is Number) ) {						throw new ODEException("Invalid style specification: " + style);					}				 this.styles[k] = style.slice();			}		}				public function clearSolutions():void {//			for (var i:int=0; i<particles.length; i++) {//				spSolutions.removeChild(particles[i]);//			}//			paricles.splice(0,particles.length);			solutionData.splice(0,solutionData.length);			initialConditions.splice(0,initialConditions.length);			spSolutions.graphics.clear();		}				public final function addSolution(x0:Vector.<Number>, t0:Number,										  t1:Number=NaN, t2:Number=NaN, timeStep:Number=NaN,										  xmin:Number=NaN, xmax:Number=NaN, 										  maxPoints:int=-1, maxChange:Number=NaN, minStep:Number=NaN,										  adaptive:Boolean =true) {			// Default values (probably reasonable)			if (isNaN(t1)) { t1 = this.tmin }			if (isNaN(t2)) { t2 = this.tmax}			if (isNaN(timeStep)) { timeStep = Math.abs(this.tmax-this.tmin)/200. }			if (isNaN(xmin)) { xmin = Math.min(this.xmin,this.ymin)-10 }			if (isNaN(xmax)) { xmax = Math.max(this.ymax,this.xmax)+10 }			if (maxPoints < 0) { maxPoints = 1000 }			if (isNaN(maxChange)) { maxChange = Math.abs(this.xmax-this.xmin)/200. }			if (isNaN(minStep)) { minStep = Math.abs(this.tmax-this.tmin)/800. }						var solutionArray:Array = solver.getSolutionAtPoints(x0, t0, t1, t2, timeStep, 																 xmin, xmax, maxPoints, maxChange, minStep, adaptive);			drawComputedSolution(x0, t0, solutionArray[0], solutionArray[1]);		}			public function drawComputedSolution(x0:Vector.<Number>, t0:Number, 											 xvalues:Array, tvalues:Array, newSolution:Boolean=true):void {			if (newSolution) {				initialConditions.push([x0.slice(0,x0.length), t0]);			}			var commands:Array = new Array();			var values:Array = new Array();			for (var j:uint=0; j<plotComponents.length; j++) {				commands[j] = new Vector.<int>();				values[j] = new Vector.<Number>();			}			if (plotType == TIMEDEPENDENT) {				for (j=0; j<plotComponents.length; j++) {					var k:uint = plotComponents[j];					commands[j].push(GraphicsPathCommand.MOVE_TO);					values[j].push(xtoPix(tvalues[0]));					values[j].push(ytoPix(xvalues[0][k]) );				}				for (j=0; j<plotComponents.length; j++) {					for (var i:int=0; i<tvalues.length; i++) {						k = plotComponents[j];						commands[j].push(GraphicsPathCommand.LINE_TO);						values[j].push(xtoPix(tvalues[i]));						values[j].push(ytoPix(xvalues[i][k]) );					}				}			} else if (plotType == PHASESPACE) {				for (j=0; j<plotComponents.length; j++) {					var k1:uint = plotComponents[j][0];					var k2:uint = plotComponents[j][1];					commands[j].push(GraphicsPathCommand.MOVE_TO);					values[j].push(xtoPix(xvalues[0][k1]));					values[j].push(ytoPix(xvalues[0][k2]) );				}				for (j=0; j<plotComponents.length; j++) {					for (i=0; i<tvalues.length; i++) {						k1 = plotComponents[j][0];						k2 = plotComponents[j][1];						commands[j].push(GraphicsPathCommand.LINE_TO);						values[j].push(xtoPix(xvalues[i][k1]));						values[j].push(ytoPix(xvalues[i][k2]));					}				}			} else {				throw new ODEException('invalid plot type');			}						if (keepSolutionData) {				solutionData.push([xvalues.slice(), tvalues.slice()]);			}						with (spSolutions.graphics) {				for (k=0; k<plotComponents.length; k++) {					style = styles[k];					lineStyle(style[0], style[1], style[2]);					drawPath(commands[k], values[k]);				}			}			//			var particle:Sprite = new Sprite();//			with (particle.graphics) {//				lineStyle(1,0x000000);//				beginFill(particleColor,0.8);//				drawCircle(0,0,particleSize);//				endFill();//			}//			particle.x = xtoPix(x0[0]);//			particle.y = ytoPix(x0[1]);//			particle.visible = showParticles_;//			particles.push(particle);//			spSolutions.addChild(particle);					}		public function recomputeSolutions():void {			var ic:Array = initialConditions.slice(0,initialConditions.length);			clearSolutions();			for (var i:int=0; i < ic.length; i++) {				var x0:Vector.<Number> = ic[i][0];				var t0:Number = ic[i][1];				addSolution(x0, t0);			}			// The time value here is arbitrary			updateArrows(t0); // The system may have changed.		}		public function getSolutionAt(k:uint, t:Number):Vector.<Number> {			if (k>solutionData.length) {				throw new ODEException("Solution at level " + k + " is not defined");			}				var times:Array = solutionData[k][0]				var data:Array = solutionData[k][1];			if (t<times[0] || t>times[times.length-1]) {				return null;			}			var dmin:Number = Infinity;			var vmin:Vector.<Number>;			for (var i:int=0; i<times.length; i++) {				var d:Number = Math.abs(t - times[i]);				if (d<dmin) {					dmin = d;					vmin = data[i];				}			}			return vmin;		}				private function onMouseClick(e:MouseEvent):void {			if (!drawSolutionOnClick) {				return;			}			var inputx:Number = xfromPix(e.localX);			var inputy:Number = yfromPix(e.localY);			var x0:Vector.<Number>			var t0:Number;			if (system.dimension==2 && plotType==PHASESPACE) {				x0 = new <Number>[inputx, inputy];				//drawSolution(x0, 0.0, tmin, tmax, (tmax-tmin)/100., 				//			 Math.min(xmin,ymin)-10, Math.max(xmax,ymax)+10, 0, (xmax-xmin)/100);				try {					addSolution(x0, 0.0);				} catch (err:ODEException) {					trace("Error plotting solution");				}			} else if (system.dimension==1 && plotType==TIMEDEPENDENT) {				x0 = new <Number>[inputy];				try {					addSolution(x0, inputx);				} catch (err:ODEException) {					trace("Error plotting solution");									}			}		}	}}