// Our only include (math is included in decide.h)
#include "decide.h"

// number of LICs which have been implemented
int LICS_USED = 15; //this should be set to 15 for final

// forward declarations
boolean LCM_computation(int i, int j);

boolean LIC_0();
boolean LIC_1();
boolean LIC_2();
boolean LIC_3();
boolean LIC_4();
boolean LIC_5();
boolean LIC_6();
boolean LIC_7();
boolean LIC_8();
boolean LIC_9();
boolean LIC_10();
boolean LIC_11();
boolean LIC_12();
boolean LIC_13();
boolean LIC_14();


boolean circle(double X1, double X2, double X3, double Y1, double Y2, double Y3, double Radius);
boolean Area(double X1, double X2, double X3, double Y1, double Y2, double Y3, double Radius);
boolean Angle(double X1, double X2, double X3, double Y1, double Y2, double Y3, double Epsilon);

void DECIDE(void) 
{
	int i,j;
	boolean temp;
	
	// set components of the CMV
	CMV[0] = LIC_0();
	CMV[1] = LIC_1();
	CMV[2] = LIC_2();
	CMV[3] = LIC_3();
	CMV[4] = LIC_4();
	CMV[5] = LIC_5();
	CMV[6] = LIC_6();
	CMV[7] = LIC_7();
	CMV[8] = LIC_8();
	CMV[9] = LIC_9();
	CMV[10] = LIC_10();
	CMV[11] = LIC_11();
	CMV[12] = LIC_12();
	CMV[13] = LIC_13();
	CMV[14] = LIC_14();

	// use CMV to compute PUM
	for (i=0; i<LICS_USED; i++) 
	{
		for (j=i+1; j<LICS_USED; j++) 
		{ // exclude diagonal with i+1
		PUM[i][j] = LCM_computation(i,j);
		PUM[j][i] = PUM[i][j]; // set other side of triangular matix
		}
	}

	// use diagonal of PUM to fill FUV
	for (i=0; i<LICS_USED; i++) 
	{
		// if PUM ignores input set FUV to 1
		if (PUM[i][i] == 0) 
		{
			FUV[i] = 1;
		} 
		else 
		{
			// otherwise check all other values in the row
			temp = 1; // gets set to 0 if the test fails
			for (j=0; j<LICS_USED; j++) 
			{
				if (PUM[i][j] == 0) 
				{
					temp = 0;
					break;
				}
			}
			FUV[i] = temp;
		}
	}

	// use FUV to set LAUNCH
	temp = 1;
	for (i=0; i<LICS_USED; i++) 
	{
		if (FUV[i] == 0) 
		{
			temp = 0;
			break;
		}
	}

	LAUNCH = temp;

}

// uses CMV to compute values for LCM
// returns 1 if the connection yields true
// returns 0 otherwise
boolean LCM_computation(int i, int j) 
{
	switch (LCM[i][j]) 
	{
		case NOTUSED:
			return 1;
		case ORR:
			return CMV[i] || CMV[j];
		case ANDD:
			return CMV[i] && CMV[j];
	}
	return 1;
}

// returns 1, (true) if at least two consecutive points are further than LENGTH1 apart
// returns 0, (false) otherwise
boolean LIC_0() 
{
	int i; 
	double xdiff, ydiff, dist;

	for (i=0; i < (NUMPOINTS-1); i++) 
	{
		xdiff = X[i] - X[i+1];
		ydiff = Y[i] - Y[i+1];
		dist = sqrt(pow(xdiff,2)+pow(ydiff,2));

		if (DOUBLECOMPARE(dist, PARAMETERS.LENGTH1) == GT) 
		{
			return 1; //return true that there are two points greater than LENGTH1
		}
	}

  return 0;//return false that there are not two points greater than LENGTH1
}

// returns 1 (true) that at least one set of three consecutive points cannot all be 
//      contained within or on a circle of radius RADIUS1
// returns 0 (false) otherwise
boolean LIC_1()
{ 	//printf("\n\nSTART LIC 1\n\n");
  	int i;
	double X1,Y1,X2,Y2,X3,Y3;
	for (i=0; i < (NUMPOINTS-2); i++) 
	{
    		//find the length between each of the three points
   		X1 = X[i];
    		Y1 = Y[i];
   		X2 = X[i+1];
    		Y2 = Y[i+1];
   		X3 = X[i+2];
    		Y3 = Y[i+2];
		if(circle(X1,X2,X3,Y1,Y2,Y3,PARAMETERS.RADIUS1) == 1)
		{
			return 1;
		}
		
	}
	//printf("\nAll points passed.\n");
	return 0;//no set of 3 points are outside of circle of RADIUS1
}



// returns 1, (true) if at least 3 consecutive points form an angle < or > than pi +/- epsilon
// returns 0, (false) otherwise
boolean LIC_2() 
{
	//printf("\n\nSTART LIC 2\n\n");
	int i;
 	double X1, X2, X3, Y1, Y2, Y3;
	for (i=0; i < (NUMPOINTS-2); i++) 
	{
   		X1 = X[i];
    		Y1 = Y[i];

   		X2 = X[i+1];
    		Y2 = Y[i+1];

   		X3 = X[i+2];
    		Y3 = Y[i+2];
		
		if(Angle(X1, X2, X3, Y1, Y2, Y3, PARAMETERS.EPSILON) == 1)
		{
			return 1;
		}
		
		
		
	}
	//printf("No angle is greater than or less than pi +/- epsilon");
	return 0;//no set of 3 points make an angle less or greater than pi+/- epsilon
}

// returns 1, (true) if at least 3 consecutive points form a triangle with angle greater than AREA1
// returns 0, (false) otherwise
boolean LIC_3() 
{
	//printf("\n\nSTART LIC 3\n\n");
	int i;
 	double X1, X2, X3, Y1, Y2, Y3;

	for (i=0; i < (NUMPOINTS-2); i++) 
	{
    		//find the length between each of the three points
   		X1= X[i];
    		Y1 = Y[i];
   		X2 = X[i+1];
    		Y2 = Y[i+1];
   		X3 = X[i+2];
    		Y3 = Y[i+2];
		
		if(Area(X1,X2,X3,Y1,Y2,Y3,PARAMETERS.AREA1) == 1)
		{
			//printf("\nA triangle has area greater than area given.\n");
			return 1; //A set of 3 points are greater triangle with area
		}
	}
	//printf("\nNo triangle has area greater than area given.\n");
	return 0;
}

//returns 1, (true) that there are a set of Q_PTS that exist in more then QUADS quadrants
//returns 0, (false) that the set of Q_PTS exist in QUADS or less.
boolean LIC_4()
{
	//printf("\n\nSTART LIC 4\n\n");
	int i,j,end,Q1, Q2, Q3, Q4;
	end = NUMPOINTS-PARAMETERS.Q_PTS;//find how far we can go with the size of the set of points we have

	for (i=0; i <= end; i++) 
	{
		//printf("i: %i,end: %i\n",i,end);
		Q1 =0;
		Q2 =0;
		Q3 =0;
		Q4 =0;
		for(j=0; j < PARAMETERS.Q_PTS; j++)//go through the set of points and map quadrants
		{
			if(
			( (DOUBLECOMPARE(X[i+j], 0) == EQ) || (DOUBLECOMPARE(X[i+j], 0) == GT) )&& 
			( (DOUBLECOMPARE(Y[i+j], 0) == EQ) || (DOUBLECOMPARE(Y[i+j], 0) == GT) )
			)
			{
				//printf("Quad 1: %f, %f\n",X[i+j],Y[i+j]);
				Q1 = 1;
				//QUAD 1
			}
			else if(
			( (DOUBLECOMPARE(X[i+j], 0) == LT) ) && 
			( (DOUBLECOMPARE(Y[i+j], 0) == EQ) || (DOUBLECOMPARE(Y[i+j], 0) == GT) )
			)
			{
				//printf("Quad 2: %f, %f\n",X[i+j],Y[i+j]);
				Q2 = 1;
				//QUAD 2
			}
			else if(
			( (DOUBLECOMPARE(X[i+j], 0) == LT) || (DOUBLECOMPARE(X[i+j], 0) == EQ) )&& 
			( (DOUBLECOMPARE(Y[i+j], 0) == LT) )
			)
			{
				//printf("Quad 3: %f, %f\n",X[i+j],Y[i+j]);
				Q3 = 1;
				//QUAD 3
			}
			else
			{
				//printf("Quad 4: %f, %f\n",X[i+j],Y[i+j]);
				Q4 = 1;
				//QUAD 4
			}
		} 
		//printf("Quads: %i\n",(Q1+Q2+Q3+Q4));
		if((Q1+Q2+Q3+Q4)>PARAMETERS.QUADS)
		{
			//printf("\nThere are more points in quadrants than quads given.\n");
			return 1;//there are sets of points that are in more quadrants than PARAMETERS.QUADS
		}
	}
	//printf("\nThere are not more points in quadrants than quads given.\n");
	return 0;
}	

// returns 1, (true) if at least two consecutive points have a difference less than 0
// returns 0, (false) otherwise
boolean LIC_5() 
{
	//printf("\n\nSTART LIC 5\n\n");
	int j; 
	double xdiff;
	for (j=1; j <= (NUMPOINTS-1); j++) 
	{
	    xdiff = X[j] - X[j-1];
	    if (DOUBLECOMPARE(xdiff, 0) == LT)
	    {
		//printf("\nThere are two points where the first X is greater than the second.\n");
		//printf("\nPoints failed: X1: %f, X2: %f\n", X[j-1],X[j]);
	    	return 1; //return true that there are two points whose difference are less than 0
	    }
	}
	//printf("\nThere are not two points where the first X is greater than the second.\n");
	return 0;//return false that there are not two points whose difference is less than 0
}

// returns 1, (true) if there is at least one set of points in N_PTS where the distance is between a point and the from first and last points is > DIST
// returns 0, (false) otherwise
boolean LIC_6() 
{
	//printf("\n\nSTART LIC 6\n\n");
	if(NUMPOINTS < 3)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}
	int i; 
	int end = NUMPOINTS-PARAMETERS.N_PTS;

	for (i=0; i <= end; i++) 
	{
		int innerpoint;
	 	double xdiff, ydiff, sideA, sideB, sideC, angleA, angleB, angleC, piHalf, height;
		piHalf = PI/2;
		for (innerpoint=0; innerpoint < (PARAMETERS.N_PTS-2); innerpoint++) 
		{
	    		//find the length between each of the three points
	   		xdiff = X[i] - X[i+PARAMETERS.N_PTS-1];
	    		ydiff = Y[i] - Y[i+PARAMETERS.N_PTS-1];
	    		sideA = sqrt(pow(xdiff,2)+pow(ydiff,2)); //distance formula

	   		xdiff = X[i+PARAMETERS.N_PTS-1] - X[i+innerpoint+1];
	    		ydiff = Y[i+PARAMETERS.N_PTS-1] - Y[i+innerpoint+1];
	    		sideB = sqrt(pow(xdiff,2)+pow(ydiff,2));

	   		xdiff = X[i] - X[i+innerpoint+1];
	    		ydiff = Y[i] - Y[i+innerpoint+1];
	    		sideC = sqrt(pow(xdiff,2)+pow(ydiff,2));
			//printf("\nSideA: %f, SideB: %f, SideC: %f\n",sideA,sideB,sideC);
			//printf("\nX1:%f, Y1:%f, X2:%f, Y2:%f, X3:%f, Y3:%f\n",X[i], Y[i], X[i+innerpoint+1] ,Y[i+innerpoint+1],X[i+PARAMETERS.N_PTS-1],Y[i+PARAMETERS.N_PTS-1]);
			if(DOUBLECOMPARE(sideA, 0) == EQ)//check that first and last point aren't 0 distance
			{
				//printf("\nEnd points the same. DIST: %f, sideB: %f\n",PARAMETERS.DIST, sideB);
				if(DOUBLECOMPARE(sideB, PARAMETERS.DIST) == GT)
				{
					//printf("\nDistance Greater than DIST.\n");
					return 1; //end points are the same and the dist between one of the points is greater than the DIST
				}
				else
				{
	
					continue;//go to next point
				}
			}
			if((DOUBLECOMPARE(sideC, 0) == EQ) ||(DOUBLECOMPARE(sideB, 0) == EQ))//check that first and last point aren't 0 distance
			{
				//printf("\nInner Point the same as end points, going to next point.\n");
				continue;//go to next point because this means the innerpoint is the same as one of the end points
			}

			//find all the angles, using the law of cosines
			double temp = ((pow(sideB,2) + pow(sideC,2) - pow(sideA,2)) / (2*sideB*sideC));
			if(DOUBLECOMPARE(temp, 1.0) == EQ)//This is necessary because of the double imprecision, you can
			{				//get 1.00000000001, that doesn't work for acos as its only -1 to 1
				angleA = acos(1.0);
			}
			else if(DOUBLECOMPARE(temp, -1.0) == EQ)
			{
				angleA = acos(-1.0);
			}
			else
			{
				angleA = acos(temp);
			}

			temp = ((pow(sideA,2) + pow(sideC,2) - pow(sideB,2)) / (2*sideA*sideC));
			if(DOUBLECOMPARE(temp, 1.0) == EQ)//This is necessary because of the double imprecision, you can
			{				//get 1.00000000001, that doesn't work for acos as its only -1 to 1
				angleB = acos(1.0);
			}
			else if(DOUBLECOMPARE(temp, -1.0) == EQ)
			{
				angleB = acos(-1.0);
			}
			else
			{
				angleB = acos(temp);
			}

			angleC = PI - (angleA+angleB); //triangles angles add up to 180 degrees
		
			//printf("\nAngleA: %f, AngleB: %f, AngleC: %f\n",angleA,angleB,angleC);
			////printf("pi is:%f, add is:%f\n",PI,(angleA+angleB+angleC));	

			//if the end angles are 90 degrees or greater then use 	the side distance
			if((DOUBLECOMPARE(angleC, piHalf) == GT) ||(DOUBLECOMPARE(angleC, piHalf) == EQ))
			{
				double newAngle = PI - angleC;
				double distToCompare = sin(newAngle)*sideB;
				//printf("\nAngle C is 90 degrees or greater.\n");
				if(DOUBLECOMPARE(distToCompare, PARAMETERS.DIST) == GT)
				{
					//printf("\nSide  greater than dist.\n");
					return 1; //distance between inner point and line made by end points is greater than DIST
				}
				else
				{
					continue;
				}
				/*if(DOUBLECOMPARE(sideB, PARAMETERS.DIST) == GT)//we use the adjacent side not the opposite of the angle as to get the closest line to the end points
				{
					//printf("\nSide B greater than dist.\n");
					return 1; //distance between inner point and line made by end points is greater than DIST
				}
				else
				{
					continue;
				}*/
			}
			else if((DOUBLECOMPARE(angleB, piHalf) == GT) ||(DOUBLECOMPARE(angleB, piHalf) == EQ))
			{
				//printf("\nAngle B is 90 degrees or greater.\n");
				double newAngle = PI - angleB;
				double distToCompare = sin(newAngle)*sideC;
				//printf("\nAngle C is 90 degrees or greater.\n");
				if(DOUBLECOMPARE(distToCompare, PARAMETERS.DIST) == GT)
				{
					//printf("\nSide  greater than dist.\n");
					return 1; //distance between inner point and line made by end points is greater than DIST
				}
				else
				{
					continue;
				}
				/*if(DOUBLECOMPARE(sideC, PARAMETERS.DIST) == GT)
				{
					//printf("\nSide C greater than dist.\n");//we use the adjacent side not the opposite of the angle
					return 1; //distance between inner point and line made by end points is greater than DIST
				}
				else
				{
					continue;
				}*/
			}			


			height = sideC*sin(angleB); //find the height of our triangle, this will be at 90 degrees, so its the sortest distance to the line made by two endpoints.
			//printf("\nHeight: %f, DIST: %f\n",height, PARAMETERS.DIST);
			if(DOUBLECOMPARE(height, PARAMETERS.DIST) == GT)
			{
				//printf("\nHeight greater than dist.\n");
				return 1; //distance between inner point and line made by end points is greater than DIST
			}

		}
	}
	//printf("\nAll points passed.\n");
	return 0;
}

// returns 1, (true) if at least 2 points separted by K_PTS are distance greater than LENGTH1
// returns 0, (false) otherwise
boolean LIC_7() 
{
	//printf("\n\nSTART LIC 7\n\n");
	if(NUMPOINTS < 3)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}

	int i; 
	double xdiff, ydiff, dist;
	int end = NUMPOINTS-1-PARAMETERS.K_PTS;

	for (i=0; i < end; i++) 
	{
		xdiff = X[i] - X[i+PARAMETERS.K_PTS+1];
		ydiff = Y[i] - Y[i+PARAMETERS.K_PTS+1];
		dist = sqrt(pow(xdiff,2)+pow(ydiff,2));
		//printf("\nX1:%f, Y1:%f, X2:%f, Y2:%f\n",X[i], Y[i], X[i+PARAMETERS.K_PTS+1] ,Y[i+PARAMETERS.K_PTS+1]);
		if (DOUBLECOMPARE(dist, PARAMETERS.LENGTH1) == GT) 
		{
			//printf("\nTwo points distance greater(%f) than LENGTH1: %f\n",dist,PARAMETERS.LENGTH1);
			return 1; //return true that there are two points greater than LENGTH1
		}
	}
	//printf("\nAll points passed.\n");
	return 0;
}

// returns 1, (true) if at least 3 points separated by A_PTS, and B_PTS can't fit in circle with RADIUS1
// returns 0, (false) otherwise
boolean LIC_8() 
{
	//printf("\n\nSTART LIC 8\n\n");
	if(NUMPOINTS < 5)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}

  	int i, end;
	double X1,Y1,X2,Y2,X3,Y3;
	end = NUMPOINTS - 2 - PARAMETERS.A_PTS - PARAMETERS.B_PTS;
	for (i=0; i < end; i++) 
	{
    		//find the length between each of the three points
   		X1 = X[i];
    		Y1 = Y[i];
   		X2 = X[i+PARAMETERS.A_PTS+1];
    		Y2 = Y[i+PARAMETERS.A_PTS+1];
   		X3 = X[i+PARAMETERS.A_PTS+PARAMETERS.B_PTS+2];
    		Y3 = Y[i+PARAMETERS.A_PTS+PARAMETERS.B_PTS+2];
		if(circle(X1,X2,X3,Y1,Y2,Y3,PARAMETERS.RADIUS1) == 1)
		{
			return 1;
		}
		
	}
	//printf("\nAll points passed.\n");
	return 0;//no set of 3 points are outside of circle of RADIUS1
}


boolean LIC_9()
{
	//printf("\n\nSTART LIC 9\n\n");
	int i;
	double X1,Y1,X2,Y2,X3,Y3;
	// check for valid NUMPOINTS
	if (NUMPOINTS < 5)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}

	// check all angles
	for (i=0; i < NUMPOINTS-(PARAMETERS.C_PTS + PARAMETERS.D_PTS)-2; i++)
	{
		X1 = X[i];
		Y1 = Y[i];

		X2 = X[i + PARAMETERS.C_PTS + 1];
		Y2 = Y[i + PARAMETERS.C_PTS + 1];

		X3 = X[i + PARAMETERS.C_PTS + PARAMETERS.D_PTS + 2];
		Y3 = Y[i + PARAMETERS.C_PTS + PARAMETERS.D_PTS + 2];
		if (Angle( X1, X2, X3, Y1, Y2, Y3, PARAMETERS.EPSILON) == 1)
		{
			return 1;
		}
	}
	// all angles failed
	//printf("\nAll points passed.\n");
	return 0;
}


boolean LIC_10()
{
	//printf("\n\nSTART LIC 10\n\n");
	int i;
	double X1,Y1,X2,Y2,X3,Y3;
	// check for valid NUMPOINTS
	if (NUMPOINTS < 5)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}

	// check all triangles
	for (i=0; i < (NUMPOINTS-(PARAMETERS.E_PTS + PARAMETERS.F_PTS)-2); i++)
	{
		X1 = X[i];
		Y1 = Y[i];

		X2 = X[i + PARAMETERS.E_PTS + 1];
		Y2 = Y[i + PARAMETERS.E_PTS + 1];

		X3 = X[i + PARAMETERS.E_PTS + PARAMETERS.F_PTS + 2];
		Y3 = Y[i + PARAMETERS.E_PTS + PARAMETERS.F_PTS + 2];
		if (Area( X1, X2, X3, Y1, Y2, Y3, PARAMETERS.AREA1) == 1)
		{
			return 1;
		}
	}
	// all triangles failed
	//printf("\nAll points passed.\n");
	return 0;
}

boolean LIC_11()
{
	//printf("\n\nSTART LIC 11\n\n");
	int i;

	// Check that NUMPOINTS is >= 3
	if (NUMPOINTS < 3)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}

	// test all pairs
	for (i=0; i < (NUMPOINTS - PARAMETERS.G_PTS - 1); i++)
	{
		if (DOUBLECOMPARE((X[i + PARAMETERS.G_PTS +  1] - X[i]), 0) == LT)
		{
			//printf("\nPoints failed on i %i: X1: %f, X2: %f\n",i, X[i],X[i + PARAMETERS.G_PTS +  1]);
			return 1;
		}
	}

	// all pairs failed
	//printf("\nAll points passed.\n");
	return 0;
}

// returns 1, (true) if at least two points separted by K_PTS are greater than LENGTH1 and two points are less than 
// length 2
// returns 0, (false) otherwise
boolean LIC_12() 
{
	//printf("\n\nSTART LIC 12\n\n");
	int i,j;
	double dist, xdiff, ydiff, dist1, dist2;
	dist1 = 0;
	dist2 = 0;

	if(NUMPOINTS < 3)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}
	dist1 = 0;
	dist2 = 0;
	j = PARAMETERS.K_PTS + 1;
	for(i=0; i < (NUMPOINTS - PARAMETERS.K_PTS - 1); i++)
	{	
		xdiff = X[i+j] - X[i];
		ydiff = Y[i+j] - Y[i];
		dist = sqrt(pow(xdiff,2)+pow(ydiff,2));
		if(DOUBLECOMPARE(dist, PARAMETERS.LENGTH1) == GT)
		{
			//printf("\nPoints found with dist: %f greater than LENGTH1: %f\n", dist, PARAMETERS.LENGTH1);
			//printf("\nPoints are:\nX1: %f, Y1: %f\nX2: %f, Y2: %f\n", X[i], Y[i],  X[i+j], Y[i+j]);
			dist1 = 1;
		}
		if(DOUBLECOMPARE(dist, PARAMETERS.LENGTH2) == LT)
		{
			//printf("\nPoints found with dist: %f less than LENGTH2: %f\n", dist, PARAMETERS.LENGTH2);
			//printf("\nPoints are:\nX1: %f, Y1: %f\nX2: %f, Y2: %f\n", X[i], Y[i],  X[i+j], Y[i+j]);
			dist2 = 1;
		}
		if((dist1 == 1) && (dist2 == 1))
		{
			//printf("\nAll points found!\n");
			return 1;
		}
	}
	//printf("\nPoints not found :(\n");
	return 0;	
}

// returns 1, (true) if at least 3 points separated by A_PTS, and B_PTS can't fit in circle with RADIUS1 and 3 points
// can fit in circle with RADIUS2
// returns 0, (false) otherwise
boolean LIC_13() 
{
	//printf("\n\nSTART LIC 13\n\n");
	int i,test1,test2;
	test1 = 0;
	test2 = 0;
	double X1,Y1,X2,Y2,X3,Y3,Radius1,Radius2;
	Radius1 = PARAMETERS.RADIUS1;
	Radius2 = PARAMETERS.RADIUS2;
	if(NUMPOINTS < 5)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}
	for(i = 0; i < (NUMPOINTS - (PARAMETERS.A_PTS + PARAMETERS.B_PTS) - 2); i++)
	{	
		X1 = X[i];
		X2 = X[i + PARAMETERS.A_PTS + 1];
		X3 = X[i + PARAMETERS.A_PTS + PARAMETERS.B_PTS + 2];
		Y1 = Y[i];
		Y2 = Y[i + PARAMETERS.A_PTS + 1];
		Y3 = Y[i + PARAMETERS.B_PTS + PARAMETERS.A_PTS + 2];
		if ((test1 != 1) && (circle(X1,X2,X3,Y1,Y2,Y3,Radius1) == 1))
		{
			//printf("\nTriangle found that wont fit in circle with radius: %f\n", Radius1);
			//printf("\nPoints are:\nX1:%f, Y1:%f\nX2:%f, Y2:%f\nX3:%f, Y3:%f\n",X1, Y1, X2 ,Y2, X3, Y3);
			test1 = 1;	
		}
		if ((test2 != 1) && (circle(X1,X2,X3,Y1,Y2,Y3,Radius2) == 0))
		{
			//printf("\nTriangle found will fit in circle with radius: %f\n",Radius2);
			//printf("\nPoints are:\nX1:%f, Y1:%f\nX2:%f, Y2:%f\nX3:%f, Y3:%f\n",X1, Y1, X2 ,Y2, X3, Y3);
			test2 = 1;
		}
		if((test1 == 1) && (test2 == 1))
		{
			//printf("\nTwo Triangles found!\n");
			return 1;
		}
	}
	//printf("\nTwo Triangles not found :(\n");
	return 0;
}

// returns 1, (true) if at least 3 points separated by E_PTS, and F_PTS have a triangle area greater than AREA1 and
// 3 points have triangle area less than AREA2
// returns 0, (false) otherwise
boolean LIC_14() 
{
	//printf("\n\nSTART LIC 14\n\n");
	int i, test3, test4;
	test3 = 0;
	test4 = 0;
	double X1, Y1, X2, Y2, X3, Y3, Area1, Area2;
	Area1 = PARAMETERS.AREA1;
	Area2 = PARAMETERS.AREA2;
	if(NUMPOINTS < 5)
	{
		//printf("\nNot enough points.\n");
		return 0;
	}
	for(i=0; i < (NUMPOINTS - (PARAMETERS.E_PTS + PARAMETERS.F_PTS) - 2); i++)
	{	
		
		X1 = X[i];
		X2 = X[i + PARAMETERS.E_PTS + 1];
		X3 = X[i + PARAMETERS.E_PTS + PARAMETERS.F_PTS + 2];
		Y1 = Y[i];
		Y2 = Y[i + PARAMETERS.E_PTS + 1];
		Y3 = Y[i + PARAMETERS.E_PTS + PARAMETERS.F_PTS + 2];
		
		if ((test3 != 1) && (Area(X1,X2,X3,Y1,Y2,Y3,Area1) == 1))
		{
			//printf("\nTriangle found with area greater than: %f\n",Area1);
			//printf("\nPoints are:\nX1:%f, Y1:%f\nX2:%f, Y2:%f\nX3:%f, Y3:%f\n",X1, Y1, X2 ,Y2, X3, Y3);
			test3 = 1;	
		}

		if ((test4 != 1) && (Area(X1,X2,X3,Y1,Y2,Y3,Area2) == 0))
		{
			//printf("\nTriangle found with area less than: %f\n",Area2);
			//printf("\nPoints are:\nX1:%f, Y1:%f\nX2:%f, Y2:%f\nX3:%f, Y3:%f\n",X1, Y1, X2 ,Y2, X3, Y3);
			test4 = 1;
		}

		if((test3 == 1) && (test4 == 1))
		{
			//printf("\nTwo Triangles found!\n");
			return 1;
		}
	}
	//printf("\nTwo Triangles not found :(\n");
	return 0;				
}

boolean circle(double X1,double X2,double X3,double Y1,double Y2,double Y3,double Radius)
{
 	double xdiff, ydiff, sideA, sideB, sideC, diameter, angleA, angleB, angleC, piHalf, circumRadius;
	circumRadius = 0;
	diameter = Radius*2; 
	piHalf = PI/2;
	//find the length between each of the three points
	xdiff = X1 - X2;
	ydiff = Y1 - Y2;
	sideA = sqrt(pow(xdiff,2)+pow(ydiff,2)); //distance formula

	xdiff = X2 - X3;
	ydiff = Y2 - Y3;
	sideB = sqrt(pow(xdiff,2)+pow(ydiff,2));

	xdiff = X1 - X3;
	ydiff = Y1 - Y3;
	sideC = sqrt(pow(xdiff,2)+pow(ydiff,2));
	//printf("\nX1:%f, Y1:%f\nX2:%f, Y2:%f\nX3:%f, Y3:%f\n",X1, Y1, X2 ,Y2, X3, Y3);
	//check if any side is greater than the diameter
	//printf("\nSideA: %f, SideB: %f, SideC: %f\n",sideA,sideB,sideC);
	if(DOUBLECOMPARE(sideA, diameter) == GT)
	{
		//printf("\nSideA is bigger than circle diameter.\n");
		return 1;//set of 3 points are outside of circle of RADIUS1
	}
	if(DOUBLECOMPARE(sideB, diameter) == GT)
	{
		//printf("\nSideB is bigger than circle diameter.\n");
		return 1;//set of 3 points are outside of circle of RADIUS1
	}
	if(DOUBLECOMPARE(sideC, diameter) == GT)
	{
		//printf("\nSideC is bigger than circle diameter.\n");
		return 1;//set of 3 points are outside of circle of RADIUS1
	}
	if((DOUBLECOMPARE(sideC, 0) == EQ)||(DOUBLECOMPARE(sideA, 0) == EQ)||(DOUBLECOMPARE(sideB, 0) == EQ))//check that we aren't going to devide by 0
	{
		//printf("\nOne or more sides is 0 length, and longest side fits in circle so continue.\n");
		return 0;
	}
	
	//find all the angles, using the law of cosines
	double temp = ((pow(sideB,2) + pow(sideC,2) - pow(sideA,2)) / (2*sideB*sideC));
	if(DOUBLECOMPARE(temp, 1.0) == EQ)//This is necessary because of the double imprecision, you can
	{				//get 1.00000000001, that doesn't work for acos as its only -1 to 1
		angleA = acos(1.0);
	}
	else if(DOUBLECOMPARE(temp, -1.0) == EQ)
	{
		angleA = acos(-1.0);
	}
	else
	{
		angleA = acos(temp);
	}

	temp = ((pow(sideA,2) + pow(sideC,2) - pow(sideB,2)) / (2*sideA*sideC));
	if(DOUBLECOMPARE(temp, 1.0) == EQ)//This is necessary because of the double imprecision, you can
	{				//get 1.00000000001, that doesn't work for acos as its only -1 to 1
		angleB = acos(1.0);
	}
	else if(DOUBLECOMPARE(temp, -1.0) == EQ)
	{
		angleB = acos(-1.0);
	}
	else
	{
		angleB = acos(temp);
	}

	angleC = PI - (angleA+angleB); //triangles angles add up to 180 degrees
	
	//printf("\nAngleA: %f, AngleB: %f, AngleC: %f\n",angleA,angleB,angleC);
	////printf("pi is:%f, add is:%f\n",PI,(angleA+angleB+angleC));		

/*
	We check if any angle is greater than or equal to 90 degrees, if so then the opposite side 
	is the hypotenuse and we center our circle on the mid point of that hypotenuse. If no angle is
	greater than or equal to 90 degrees then we find the radius of the circumcircle and check if its
	bigger than the radius we want.
*/
	if((DOUBLECOMPARE(angleC, piHalf) == GT) ||(DOUBLECOMPARE(angleC, piHalf) == EQ))
	{
		//printf("\nAngle C is 90 degrees or greater.\n");
		return 0;//if obtuse or right then we can fit as we already checked that the hypotenuse is less than the diameter.
	}
	else if((DOUBLECOMPARE(angleB, piHalf) == GT) ||(DOUBLECOMPARE(angleB, piHalf) == EQ))
	{
		//printf("\nAngle B is 90 degrees or greater.\n");
		return 0;
	}
	else if((DOUBLECOMPARE(angleA, piHalf) == GT) ||(DOUBLECOMPARE(angleA, piHalf) == EQ))
	{
		//printf("\nAngle A is 90 degrees or greater.\n");
		return 0;
	}
	else //acute triangles can use the circumcircles radius and if its bigger than RADIUS1 return 1
	{

		circumRadius = (sideA*sideB*sideC)/sqrt((sideA+sideB+sideC)*(sideC+sideB-sideA)*(sideA+sideC-sideB)*(sideA+sideB-sideC));
		//printf("\ncircumR is: %f\n",circumRadius);
		if(DOUBLECOMPARE(circumRadius, Radius) == GT)
		{
			//printf("\nAcute triangle wont fit\n");
			return 1;//one of the three points doesn't exist in the circle	
		}
		//printf("\nAcute triangle will fit\n");
	}
	/*Removed Section*/
	//printf("\nAll points passed.\n");
	return 0;//no set of 3 points are outside of circle of RADIUS1
}

boolean Area(double X1,double X2,double X3,double Y1,double Y2,double Y3,double Area)
{
	//printf("\nX1:%f, Y1:%f\nX2:%f, Y2:%f\nX3:%f, Y3:%f\n",X1, Y1, X2 ,Y2, X3, Y3);
 	double xdiff, ydiff, sideA, sideB, sideC, p, area;
	//find the length between each of the three points
	xdiff = X1 - X2;
	ydiff = Y1 - Y2;
	sideA = sqrt(pow(xdiff,2)+pow(ydiff,2));//distance formula

	xdiff = X2 - X3;
	ydiff = Y2 - Y3;
	sideB = sqrt(pow(xdiff,2)+pow(ydiff,2));

	xdiff = X1 - X3;
	ydiff = Y1 - Y3;
	sideC = sqrt(pow(xdiff,2)+pow(ydiff,2));
	
	//Find the area of the triangle using heron's formula.
	p = ((sideA+sideB+sideC)/2);
	double temp = (p*(p-sideA)*(p-sideB)*(p-sideC));

	//printf("\nTEMP IS: %f\n",temp);
	if((DOUBLECOMPARE(temp, -0) == EQ))//This is to remove -0 from messing up sqrt
	{
		temp = 0;
	}
	area = sqrt(temp);
	//printf("\n\nSTART OF AREA\n");	
	//printf("\nSideA: %f, SideB: %f, SideC: %f\n",sideA,sideB,sideC);		
	//printf("\nP is: %f, Area is: %f, Given Area: %f\n",p,area, Area);

	if((DOUBLECOMPARE(p, sideB) == EQ)||(DOUBLECOMPARE(p, sideA) == EQ)||(DOUBLECOMPARE(p, sideC) == EQ))//This is a flat line since have the perimeter is the same size as one side
	{
		double area4 = (X1*(Y2-Y3))+(X2*(Y3-Y1))+(X3*(Y1-Y2));
		area4 = area4*area4;
		area4 = sqrt(area4);
		area4 = .5*area4;
		if(DOUBLECOMPARE(area4, Area) == GT)
		{
			return 1; //distance between inner point and line made by end points is greater than DIST
		}
		return 0;
	}

	if(DOUBLECOMPARE(area, Area) == GT)
	{
		//printf("\nA triangle has area greater than area given.\n");
		return 1; //A set of 3 points are greater triangle with area
	}
	return 0;
} 

boolean Angle(double X1,double X2,double X3,double Y1,double Y2,double Y3, double Epsilon)
{
	//printf("\nX1:%f, Y1:%f\nX2:%f, Y2:%f\nX3:%f, Y3:%f\n",X1, Y1, X2 ,Y2, X3, Y3);
	double xdiff, ydiff, sideA, sideB, sideC, angleC;
	//find the length between each of the three points
	xdiff = X1 - X2;
	ydiff = Y1 - Y2;
	sideA = sqrt(pow(xdiff,2)+pow(ydiff,2));//distance formula

	xdiff = X2 - X3;
	ydiff = Y2 - Y3;
	sideB = sqrt(pow(xdiff,2)+pow(ydiff,2));

	xdiff = X1 - X3;
	ydiff = Y1 - Y3;
	sideC = sqrt(pow(xdiff,2)+pow(ydiff,2));
	
	//printf("\nSideA: %f, SideB: %f, SideC: %f\n",sideA,sideB,sideC);
	

	if((DOUBLECOMPARE(sideA, 0) == EQ)||(DOUBLECOMPARE(sideB, 0) == EQ))
	{
		//side a or b is 0 length, side c can be 0 length as its opposite the vertex
		//printf("\nA point is same as the vertex.\n");
		return 0;
	}
	

	//find the vertex using the law of cosines
	double temp = ((pow(sideB,2) + pow(sideA,2) - pow(sideC,2)) / (2*sideB*sideA));
	if(DOUBLECOMPARE(temp, 1.0) == EQ)//This is necessary because of the double imprecision, you can
	{				//get 1.00000000001, that doesn't work for acos as its only -1 to 1
		angleC = acos(1.0);
	}
	else if(DOUBLECOMPARE(temp, -1.0) == EQ)
	{
		angleC = acos(-1.0);
	}
	else
	{
		angleC = acos(temp);
	}		
		
	
	//printf("\nAngleC: %f, Angle given: %f, %f\n",angleC,(PI + Epsilon),(PI- Epsilon));
	
	if((DOUBLECOMPARE(angleC, (PI + Epsilon)) == GT) ||(DOUBLECOMPARE(angleC, (PI - Epsilon)) == LT))
	{
		//printf("An angle is greater than or less than pi +/- epsilon");
		return 1; //A set of 3 points are greater or less than pi +/- epsilon
	}
	//printf("\nSideA: %f, SideB: %f, SideC: %f\n",sideA,sideB,sideC);
	return 0;
}
