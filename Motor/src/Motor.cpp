/*  Simulateur de physique : orbites planétaires
 *      Auteur : Maxence Lemoine
 *      Date : 28/01/2025
 *
 */

#include "pch.h"
//#include "pgl.h"

int main(int argc, char** argv) {

	Application app(500, 500, "TIPE");

	while (!app.ShouldClose())
	{
		app.Run();
	}

	return 0;
}