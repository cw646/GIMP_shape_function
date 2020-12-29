#include <iostream>

#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <iomanip>      // std::setprecision

using Eigen::MatrixXd;
int main() {

    Eigen::Matrix<double, 2, 1> xi;
    Eigen::Matrix<double, 2, 1> particle_size;

    //! Dimension
    const unsigned Dim = 2;
    //! Nodes in GIMP function
    // 16 for 2D 4 for 1D
    const unsigned Nfunctions = 4;

    //! To store shape functions
    Eigen::Matrix<double, 16, 1> shapefn;

    //! To store grad shape functions
    Eigen::Matrix<double, 16, 2> grad_shapefn;

    //! length of element in local coordinate
    // Keep this value constant
    const double element_length = 2.;

    // Local coordinates of 2D GIMP cells
    /*const Eigen::Matrix<double, Nfunctions, Dim> local_nodes =
            (Eigen::Matrix<double, Nfunctions, Dim>() << -1., -1.,

                    1., -1.,
                    1.,  1.,
                    -1.,  1.,
                    -3., -3.,
                    -1., -3.,
                    1., -3.,
                    3., -3.,
                    3., -1.,
                    3.,  1.,
                    3.,  3.,
                    1.,  3.,
                    -1.,  3.,
                    -3.,  3.,
                    -3.,  1.,
                    -3., -1.).finished();*/

    // Local coordinates of 1D GIMP cells
    const Eigen::Matrix<double, Nfunctions, Dim> local_nodes =
            (Eigen::Matrix<double, Nfunctions, Dim>() << 

		-3.,0.,
		-1.,0.,
		 1.,0.,
		 3.,0.
		).finished();


    // particle location
    xi << -4.0, 0.0;
    // particle size
    particle_size << 1.1, 1.1;

    double interval = 0.05;

    std::vector<double> function_store;
    std::vector<double> gfunction_store;

    
    for (unsigned k = 0; k < 120; ++k) {
	
	xi(0) += interval;

	//Function loop
        for (unsigned n = 0; n < Nfunctions; ++n) {
            Eigen::Matrix<double, 2, 1> sni;
            Eigen::Matrix<double, 2, 1> dni;

		// GIMP conditional statement loop
            for (unsigned i = 0; i < Dim; ++i) {
		//length of particle
                double lp = particle_size(i) * 0.5;
		// active node
                double ni = local_nodes(n, i);
		// local particle  - local node
                double npni = xi(i) - ni;  
                //! Conditional shape function statement
                // see: Pruijn, N.S., 2016. Eq(4.30)
                if (npni <= (-element_length - lp)) {
                    sni(i) = 0.;
                    dni(i) = 0.;
                } else if ((-element_length - lp) < npni &&
                           npni <= (-element_length + lp)) {

                    sni(i) = std::pow(element_length + lp + npni, 2.) /
                             (4. * (element_length * lp));
                    dni(i) = (element_length + lp + npni) / (2. * element_length * lp);
                } else if ((-element_length + lp) < npni && npni <= -lp) {
                    sni(i) = 1. + (npni / element_length);
                    dni(i) = 1. / element_length;
                } else if (-lp < npni && npni <= lp) {
                    sni(i) =
                            1. - (((npni * npni) + (lp * lp)) / (2. * element_length * lp));
                    dni(i) = -(npni / (element_length * lp));
                } else if (lp < npni && npni <= (element_length - lp)) {
                    sni(i) = 1. - (npni / element_length);
                    dni(i) = -(1. / element_length);
                } else if ((element_length - lp) < npni &&
                           npni <= (element_length + lp)) {
                    sni(i) = std::pow(element_length + lp - npni, 2.) /
                             (4. * element_length * lp);
                    dni(i) = -((element_length + lp - npni) / (2. * element_length * lp));
                } else if ((element_length + lp) < npni) {
                    sni(i) = 0.;
                    dni(i) = 0.;
                } else {
                    throw std::runtime_error(
                            "GIMP grad shapefn: Point location outside area of influence");
                }
            }
            // 2D Shape
            //shapefn(n) = sni(0) * sni(1);

            // 1D Shape value @ node n
	    shapefn(n) = sni(0);
		
	    // store 1D shape function value at specific node (n)
	    if (n == 1 )
	    function_store.push_back(shapefn(n));
           

            // 2D Grad
            //grad_shapefn(n, 0) = dni(0) * sni(1);
            //grad_shapefn(n, 1) = dni(1) * sni(0);
            
            // 1D Grad
            grad_shapefn(n, 0) = dni(0);
	    
	    // store 1D grad function value at specific node (n)
	    if (n == 1 )
	    gfunction_store.push_back(grad_shapefn(n,0));
        }
    }
    std::cout << "function store size: " << function_store.size() << '\n';
    for (unsigned n = 0; n < Nfunctions; ++n) {
        std::cout << n << " "
                  << " s(" << local_nodes(n, 0) << " , " << local_nodes(n, 1)
                  << ") : " << shapefn(n) << '\n';
        //std::cout << n << " " << std::setprecision(9) << grad_shapefn(n,1) << '\n';
    }

    //! Output file
    std::string xvalfilename = "snvals.txt";
    std::fstream xvalfile;
    xvalfile.open(xvalfilename, std::ios::out);

    if (xvalfile.is_open()) {
        //! Write
	//function_store for shape function, gfunction_store for gradient
        for (auto const& xvals : gfunction_store) {
            xvalfile << xvals << '\n';
        }
        xvalfile.close();
    }
}
