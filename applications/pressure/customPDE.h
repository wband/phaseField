#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};
    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC, const std::vector<double> &data);

    // Function to set the non-uniform Dirichlet boundary conditions (in ICs_and_BCs.h)
    void setNonUniformDirichletBCs(const dealii::Point<dim> &p, const unsigned int index, const unsigned int direction, const double time, double & scalar_BC, dealii::Vector<double> & vector_BC);

    private:
    #include "../../include/typeDefs.h"

    const userInputParameters<dim> userInputs;

    // Function to set the RHS of the governing equations for explicit time dependent equations (in equations.h)
    void explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the RHS of the governing equations for all other equations (in equations.h)
    void nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set the LHS of the governing equations (in equations.h)
    void equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                     dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

    // Function to set postprocessing expressions (in postprocess.h)
    #ifdef POSTPROCESS_FILE_EXISTS
    void postProcessedFields(const variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
                    variableContainer<dim,degree,dealii::VectorizedArray<double> > & pp_variable_list,
                    const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;
    #endif


	// ================================================================
	// Methods specific to this subclass
	// ================================================================
    const static unsigned int CIJ_tensor_size =2*dim-1+dim/3;
    dealii::Tensor<2,CIJ_tensor_size> CIJ_base0 =
            userInputs.get_model_constant_elasticity_tensor("CIJ_base0");
    // const double coupling = userInputs.get_model_constant_double("coupling");
    // const double kwell = userInputs.get_model_constant_double("kwell");
    // const double mu_bar = userInputs.get_model_constant_double("mu_bar");
    // const double L_int = userInputs.get_model_constant_double("L_int");
	// ================================================================
    dealii::Tensor<4,dim>
            VoigtToFourthRank(const dealii::Tensor<2,CIJ_tensor_size> & voigtTensor) {
            
        dealii::Tensor<4,dim> fourthRankTensor;
        if (dim == 3){
        
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                for (int k = 0; k < dim; ++k) {
                    for (int l = 0; l < dim; ++l) {
                        int indV1 = ((i == j) ? i : 6-i-j);
                        int indV2 = ((k == l) ? k : 6-k-l);
                        fourthRankTensor[i][j][k][l] = voigtTensor[indV1][indV2];
                    }
                }
            }
        }
        } else {
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                for (int k = 0; k < dim; ++k) {
                    for (int l = 0; l < dim; ++l) {
                        int indV1 = ((i == j) ? i : 3-i-j);
                        int indV2 = ((k == l) ? k : 3-k-l);
                        fourthRankTensor[i][j][k][l] = voigtTensor[indV1][indV2];
                    }
                }
            }
        }
        }
        return fourthRankTensor;
    }
    
    
    // ================================================================
	// Model constants specific to this subclass
	// ================================================================
    const dealii::Tensor<4,dim> Cijkl0 = VoigtToFourthRank(CIJ_base0);
    const double pressure = userInputs.get_model_constant_double("pressure");
};
