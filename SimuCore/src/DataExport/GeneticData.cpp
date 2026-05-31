#include <pch.h>
#include <DataExport/GeneticData.h>

std::string GeneticData::string() const {
	std::ostringstream oss;
	oss << std::setprecision(std::numeric_limits<double>::max_digits10);
	/* On veut envoyer :
	*	- numero generation
	*	- taille population
	*	- nombre genes
	*	- meilleur genes
	*	- meilleur score
	*/

	throw std::runtime_error("la fonction string de GeneticData n'est pas encore implémentee");

	return oss.str();
}