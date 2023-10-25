// =================================================================================
// NUCLEATION FUNCTIONS
// =================================================================================

// =================================================================================
// Nucleation probability
// =================================================================================
template <int dim, int degree>
double customPDE<dim,degree>::getNucleationProbability(variableValueContainer variable_value, double dV, dealii::Point<dim> p, unsigned int variable_index) const
{
	//Supersaturation factor
    double k1 = 500;
    double k2 = 5;
    double tau = 5;
    double ssf;
    if (dim ==2) ssf=variable_value(0)-cmin[0][0];
    if (dim ==3) ssf=(variable_value(0)-cmin[0][0])*(variable_value(0)-cmin[0][0]);
	// Calculate the nucleation rate
	double J=k1*exp(-k2/(std::max(ssf,1.0e-6)))*exp(-tau/(this->currentTime));
	double retProb=1.0-exp(-J*userInputs.dtValue*((double)userInputs.steps_between_nucleation_attempts)*dV);
    return retProb;
}
