package ode {
	import flash.display.Shape;
	import flash.display.DisplayObjectContainer;
	
	/*
	 Class: Arrow
	*/
	internal class Arrow extends Shape{
		public static const REF_BASE:uint     = 0;
		public static const REF_MIDDLE:uint   = 1;
		public static const REF_MIDSHAFT:uint = 2;

		private var shaftWidth:uint;
		private var pointWidth:uint;
		private var pointHeight:uint;
		private var color_:uint;
		private var thickness:uint;
		private var reference:uint;
		private var isArrow_:Boolean = true;
		
		public function Arrow(shaftWidth:Number, pointWidth:Number, pointHeight:Number, 
							  thickness:Number=1, color:uint=0x000000, 
							  reference:uint=REF_BASE) {
			this.shaftWidth = shaftWidth;
			this.pointWidth = pointWidth;
			this.pointHeight = pointHeight;
			this.color_ = color;
			this.thickness = thickness;
			this.reference = reference;
			drawArrow();
		}
		
		private function drawArrow():void {
			if (isArrow_) {
				with (graphics) {
					clear();
					lineStyle(thickness, color_, 1, true);
					if (reference==REF_BASE) {
						moveTo(0,0);
						lineTo(shaftWidth,0);
						beginFill(color_);
						moveTo(shaftWidth,-pointHeight/2);
						lineTo(shaftWidth,pointHeight/2);
						lineTo(shaftWidth+pointWidth,0);
						lineTo(shaftWidth,-pointHeight/2);
						endFill();
					} else if (reference==REF_MIDDLE) {
						var aux1:Number = (shaftWidth+pointWidth)/2;
						var aux2:Number = (shaftWidth-pointWidth)/2
						moveTo(-aux1,0);
						lineTo(aux2,0);
						beginFill(color_);
						moveTo(aux2,-pointHeight/2);
						lineTo(aux2,pointHeight/2);
						lineTo(aux1,0);
						lineTo(aux2,-pointHeight/2);
						endFill();
					} else if (reference==REF_MIDSHAFT) {
						var aux:Number = shaftWidth/2;
						moveTo(-aux,0);
						lineTo(aux,0);
						beginFill(color_);
						moveTo(aux,0);
						lineTo(aux,-pointHeight/2);
						lineTo(aux,pointHeight/2);
						lineTo(aux+pointWidth,0);
						lineTo(aux,-pointHeight/2);
						endFill();
					} 
				} 
			} else {
				with (graphics) {
					clear();
					moveTo(0,0);
					lineStyle(thickness,color_,1,true)
					beginFill(color_);
					drawCircle(0,0,shaftWidth/2);
					endFill();
				}
			}
		}
						
		
		public function set color(cval:uint):void {
			color_ = cval;
			drawArrow();
		}
		
		public function get color():uint {
			return color_;
		}
		
		public function set isArrow(b:Boolean):void {
			isArrow_ = b;
			drawArrow();
		}
		
		public function get isArrow():Boolean {
			return isArrow_;
		}
		
		public function setShape(shaftWidth:Number, pointWidth:Number, height:Number, 
								 thickness:Number, reference:uint):void {
			this.shaftWidth = shaftWidth;
			this.pointWidth = pointWidth;
			this.pointHeight = height;
			this.thickness = thickness;
			this.reference = reference;
			drawArrow();
		}
	}
}