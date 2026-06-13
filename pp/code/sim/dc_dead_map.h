bool DCdeadArea(float alpha, float phi, float board, float zed){
  
        
       if(phi > 1.5) // DCE
	{     
              if(!(alpha < 0.359 * board - 12.558 || alpha > 0.376 * board - 11.933)) return true;
              if(!(alpha > -0.779 * board + 31.108 || alpha < -0.488 * board + 18.641)) return true;
              if( alpha > 0 && !(alpha < 0.675*board - 0.725) ) return true;
	      if( alpha <0 && !(alpha > -0.124*board - 0.109) ) return true;
	      if( alpha > 0 && !(alpha < -0.136*board + 10.849) ) return true;
	      if( alpha <0 && !(alpha > 0.711*board - 55.961) ) return true;
	      
	      if(zed >= 0){   if(!(alpha < -0.417 * board + 22.752 || alpha > -0.568 * board + 31.460)) return true; }
              else {  if(!(alpha < 0.272 * board - 6.907 || alpha > 0.338 * board - 8.308)) return true; }
              
		}
		
	else // DCW
		{
			
      		        if(!(alpha > 0.70 * board - 0.70 * 33.0 || alpha < 0.70 * board - 0.70 * 34.5)) return true;
      		        if(!(alpha < -0.36 * board + 0.36 * 37.5 || alpha > -0.36 * board + 0.36 * 43.5)) return true;
      		        if(!(alpha > 0.17 * board - 0.17 * 38.5 || alpha < 0.17 * board - 0.17 * 43.5)) return true;
      		        if(!(alpha > 0.45 * board - 0.45 * 56.0 || alpha < 0.45 * board - 0.45 * 57.5)) return true;
			if( alpha <0 && !(alpha > -0.231*board + 0.167) ) return true;
      		        if( alpha >0 && !(alpha < 0.392*board - 0.444) ) return true;
			if( alpha >0 && !(alpha < -0.277*board + 21.917)) return true;
			if( alpha <0 && !(alpha > 0.135*board - 10.752) ) return true;
			
			if(zed >= 0) 
			{
			        if(!(alpha > 0.45 * board - 0.45 * 25.0 || alpha < 0.45 * board - 0.45 * 26.3)) return true;
				if(!(alpha < -0.36 * board + 0.36 * 20.5 || alpha > -0.36 * board + 0.36 * 22.0)) return true;
        		        if(!(alpha > 0.45 * board - 0.45 * 25.0 || alpha < 0.45 * board - 0.45 * 26.5)) return true;
				if(!(alpha < -0.36 * board + 0.36 * 57.0 || alpha > -0.36 * board + 0.36 * 58.5)) return true;
        		        if(!(alpha < -0.36 * board + 0.36 * 55.0 || alpha > -0.36 * board + 0.36 * 56.0)) return true;
        		        if(!(alpha < -0.36 * board + 0.36 * 35.5 || alpha > -0.36 * board + 0.36 * 39.5)) return true;
				if(!(alpha < 0.159 * board - 11.217 || alpha > 0.163 * board - 11.389)) return true;
				if(!(alpha < -0.273 * board + 4.644 || alpha > -0.311 * board + 5.596)) return true;
				if(!(alpha > 0.497 * board -8.468 || alpha < 0.447 * board - 8.103)) return true;
			}
			else 
			{
			        if(!(alpha > 0.643 * board - 5.979 || alpha < 0.467 * board - 4.676)) return true;
        		        if(!(alpha < -0.36 * board + 0.36 * 20.5 || alpha > -0.36 * board + 0.36 * 22.0)) return true;
        		        if(!(alpha > 0.45 * board - 0.45 * 25.0 || alpha < 0.45 * board - 0.45 * 26.5)) return true;
				
        		     
			}
		}// areas rejection:end
	
	// If ran till here, the region was not dead => return false
	return false;
}
