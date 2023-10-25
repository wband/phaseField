#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
class customPDE: public MatrixFreePDE<dim,degree>
{
public:
	customPDE(userInputParameters<dim> _userInputs): MatrixFreePDE<dim,degree>(_userInputs) , userInputs(_userInputs) {};
    // Function to set the initial conditions (in ICs_and_BCs.h)
    void setInitialCondition(const dealii::Point<dim> &p, const unsigned int index, double & scalar_IC, dealii::Vector<double> & vector_IC);

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
    dealii::Tensor<2,CIJ_tensor_size> CIJ_base =
            userInputs.get_model_constant_elasticity_tensor("CIJ_base");
    const double coupling = userInputs.get_model_constant_double("coupling");
    const double kwell = userInputs.get_model_constant_double("kwell");
    const double mu_bar = userInputs.get_model_constant_double("mu_bar");
    const double L_int = userInputs.get_model_constant_double("L_int");
	// ================================================================
    dealii::Tensor<4,dim>
            VoigtToFourthRank(const dealii::Tensor<2,CIJ_tensor_size> & voigtTensor) {

        dealii::Tensor<4,dim> fourthRankTensor;
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
        return fourthRankTensor;
    }
    
    dealii::Tensor<2, dim> computeRotationMatrix(double roll_deg, double pitch_deg, double yaw_deg) {
        // Convert Euler angles from degrees to radians
        const double roll_rad = roll_deg * M_PI / 180.0;
        const double pitch_rad = pitch_deg * M_PI / 180.0;
        const double yaw_rad = yaw_deg * M_PI / 180.0;

        // Compute trigonometric values
        const double cr = cos(roll_rad);
        const double sr = sin(roll_rad);
        const double cp = cos(pitch_rad);
        const double sp = sin(pitch_rad);
        const double cy = cos(yaw_rad);
        const double sy = sin(yaw_rad);

        // Compute the rotation matrix
        dealii::Tensor<2, dim> rotation_matrix;
    
        rotation_matrix[0][0] = cp * cy;
        rotation_matrix[0][1] = cp * sy;
        rotation_matrix[0][2] = -sp;
        rotation_matrix[1][0] = sr * sp * cy - cr * sy;
        rotation_matrix[1][1] = sr * sp * sy + cr * cy;
        rotation_matrix[1][2] = sr * cp;
        rotation_matrix[2][0] = cr * sp * cy + sr * sy;
        rotation_matrix[2][1] = cr * sp * sy - sr * cy;
        rotation_matrix[2][2] = cr * cp;

        return rotation_matrix;
    }
    
    Tensor<4, dim> rotateStiffnessTensor(
            const Tensor<4, dim> & stiffness_tensor,
            const Tensor<2, dim> & rotation_matrix) {

    Tensor<4, dim> rotated_stiffness_tensor;
    for (unsigned int i = 0; i < dim; ++i) {
        for (unsigned int j = i; j < dim; ++j) {
            for (unsigned int k = 0; k < dim; ++k) {
                for (unsigned int l = k; l < dim; ++l) {
                    rotated_stiffness_tensor[i][j][k][l] = 0.0;
                    for (unsigned int p = 0; p < dim; ++p) {
                        for (unsigned int q = 0; q < dim; ++q) {
                            for (unsigned int r = 0; r < dim; ++r) {
                                for (unsigned int s = 0; s < dim; ++s) {
                                    rotated_stiffness_tensor[i][j][k][l] +=
                                        rotation_matrix[i][p] * rotation_matrix[j][q] *
                                        rotation_matrix[k][r] * rotation_matrix[l][s] *
                                        stiffness_tensor[p][q][r][s];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return rotated_stiffness_tensor;
    }
    
    // ================================================================
	// Model constants specific to this subclass
	// ================================================================
    const dealii::Tensor<4,dim> Cijkl0 = VoigtToFourthRank(CIJ_base);
    const dealii::Tensor<2,dim> rot1 = computeRotationMatrix(45,45,0);
    const dealii::Tensor<4,dim> Cijkl1 = rotateStiffnessTensor(Cijkl0, rot1);
};
