﻿/*  Simulateur de physique : orbites planétaires
 *      Auteur : Maxence Lemoine
 *      Date : 28/01/2025
 *
 */

#include "pch.h"
//#include "pgl.h"

int main(int argc, char** argv) {

	pgl::Application app(640, 480, "TIPE");

	while (!app.ShouldClose())
	{
		app.Run();
	}

	return 0;
}