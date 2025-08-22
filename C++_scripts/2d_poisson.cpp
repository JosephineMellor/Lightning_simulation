#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <array>
#include <filesystem>
#include <tuple>
#include <Eigen/Sparse>

const int I0 = 106405;
const int alpha = 22708;
const int beta = 1294530;
const double PI = 3.141592653589793;
const double R0 = 2e-3; //2mm in meters
const double T0 = 298;

//defining some useful types
typedef std::array<double , 500> data_vec;
typedef std::vector<double> big_vector;
typedef std::vector<std::vector<double>> data_table;
typedef std::array<double,4> array;
typedef std::vector<std::vector<std::array<double , 4>>> big_array;
typedef std::array<double , 19> species_vec;
typedef std::array<data_table, 19> species_tables;



//function to read in the data from the tabulated equation of state
std::tuple<data_vec , data_vec , data_table , data_table , data_table , data_table , data_table, species_vec , species_tables> Plasma19(){
    
    //load the file
    std::ifstream file("Plasma19_EoS.txt"); 

    //set up the lines, starting row, and data arrays
    std::string line;
    int row = 0;
    data_vec densities;
    data_vec pressures;
    species_vec heats_of_formation;
    species_tables mass_fractions;

    //we want to save sound speeds, energies, temperatures, electrical conductivity and thermal conductivity in tables
    data_table sound_speeds(500, std::vector<double>(500));
    data_table energies(500, std::vector<double>(500));
    data_table temperatures(500, std::vector<double>(500));
    data_table electrical_conductivity(500, std::vector<double>(500));
    data_table thermal_conductivity(500, std::vector<double>(500));
    
    //now loop over the data file row by row, we only want the fifth and sixth row for densities and pressures
    while(std::getline(file , line)){
        row++;

        if( row ==6){
            std::istringstream iss(line); //extract the line

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            } //put it temporarily in an array

            //then save the heats of formation

            for (int i = 0; i < 19; ++i) {
                if (i < (int)row_data.size()) {
                    heats_of_formation[i] = row_data[i];
                } else {
                    std::cerr << "Insufficient data in row 6 for heats_of_formation\n";
                    exit(1);
                }
            }
        }

        if( row ==10 || row ==12 ){
            std::istringstream iss(line); //extract the line

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            } //put it temporarily in an array

            //then save the densities and pressures used for indexing the tables

            for(int i=0; i<500; ++i){
                if(row ==10){
                    densities[i] = row_data[i];
                }
                if(row ==12){
                    pressures[i] = row_data[i];
                }
            }
        }

        //sound speed
        if( row >= 14 && row< 514){
            std::istringstream iss(line); 

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            sound_speeds[row-14] = std::move(row_data);
        }

        //internal energy
        if( row >= 515 && row< 1015){
            std::istringstream iss(line); 

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }
            if(row_data.size() != 500) {
                std::cerr << "Row " << row << " for energy has " << row_data.size() << " values, expected 500.\n";
                exit(1);
            }

            energies[row-515] = std::move(row_data);
        }

        //temperatures
        if( row >= 1016 && row< 1516){
            std::istringstream iss(line);

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            temperatures[row-1016] = std::move(row_data);
        }

        //electrical conductivity
        if( row >= 1517 && row< 2017){
            std::istringstream iss(line); 

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            electrical_conductivity[row-1517] = std::move(row_data);
        }

        //mass fractions
        for (int i = 0; i < 19; i++) {
            int min_row = 501 * (i + 4) + 14;
            int max_row = min_row + 500;
            mass_fractions[i].resize(500, std::vector<double>(500));

            if (row >= min_row && row < max_row) {
                std::istringstream iss(line);

                big_vector row_data;
                double val;
                while (iss >> val) {
                    row_data.push_back(val);
                }

                mass_fractions[i][row - min_row] = std::move(row_data);
            }
        }


        //thermal conductivity
        if (row >= 12539 && row < 13039) {
            std::istringstream iss(line);

            big_vector row_data;
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            if(row_data.size() != 500){
                std::cerr << "Row " << row << " for thermal conductivity has " << row_data.size() << " values, expected 500.\n";
                exit(1);
            }

            thermal_conductivity[row - 12539] = std::move(row_data);
        }

    }



    return{densities , pressures , sound_speeds , energies , temperatures , electrical_conductivity , thermal_conductivity, heats_of_formation, mass_fractions};
}


//import our data tables and arrays as global variables

auto EoS_data = Plasma19(); // the tuple
const auto& densities = std::get<0>(EoS_data);
const auto& pressures = std::get<1>(EoS_data);
const auto& sound_speeds = std::get<2>(EoS_data);
const auto& energies = std::get<3>(EoS_data);
const auto& temperatures = std::get<4>(EoS_data);
const auto& electrical_conductivity = std::get<5>(EoS_data);
const auto& thermal_conductivity = std::get<6>(EoS_data);
const auto& heats_of_formation = std::get<7>(EoS_data);
const auto& mass_fractions = std::get<8>(EoS_data);

//function to perform bilinear interpolation on a simplified grid (0,0) , (0,1) , (1,0) , (1,1)
double BilinearInterpolation(double f_NE, double f_SE, double f_SW, double f_NW, double x, double y){
    double f_x_y1 = (1-x)*f_SW + x*f_SE;
    double f_x_y2 = (1-x)*f_NW + x*f_NE;
    double f_x_y = (1-y)*f_x_y1 + y*f_x_y2;

    return f_x_y;
}

//function to 'reverse' bilinear interpolation on a simplified grid (0,0) , (0,1) , (1,0) , (1,1)
double ReverseBilinearInterpolation(double f_NE, double f_SE, double f_SW, double f_NW, double x, double f_x_y){
    double f_x_y1 = (1-x)*f_SW + x*f_SE;
    double f_x_y2 = (1-x)*f_NW + x*f_NE;
    double y = (f_x_y - f_x_y1) / (f_x_y2 - f_x_y1);
    if(f_x_y2 == f_x_y1){y = 1;}

    return y;
}

//functions to convert between primitive to conseravative variables
array PrimativeToConservative(array prim){
    array consv;

    //start with the variables that dont require the equation of state
    consv[0] = prim[0]; //rho
    consv[1] = prim[0] * prim[1]; //mom_x
    consv[2] = prim[0] * prim[2]; //mom_y

    //find index for density with binary search
    int point1 =0;
    int point2 =499;
    int half;
    int density_index;
    double rho = prim[0];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(rho > densities[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up the linear interpolation
    double density_top = densities[point2];
    int density_top_i = point2;
    double density_bottom = densities[point1];
    int density_bottom_i = point1;
    double density_ratio = (rho - density_bottom) / (density_top - density_bottom);
    if(density_bottom == density_top){density_ratio =1;}
    

    //find index for pressure also with binary search
    point1 =0;
    point2 =499;
    half;
    int pressure_index;
    double p = prim[3];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(p > pressures[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up linear interpolation
    double pressure_top = pressures[point2];
    int pressure_top_i = point2;
    double pressure_bottom = pressures[point1];
    int pressure_bottom_i = point1;
    double pressure_ratio = (p - pressure_bottom) / (pressure_top - pressure_bottom);
    if(pressure_bottom == pressure_top){pressure_ratio =1;}

    //then apply to the energies table
    double energies_NE = energies[pressure_top_i][density_top_i];
    double energies_SE = energies[pressure_bottom_i][density_top_i];
    double energies_SW = energies[pressure_bottom_i][density_bottom_i];
    double energies_NW = energies[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double energy = BilinearInterpolation(energies_NE , energies_SE , energies_SW , energies_NW , density_ratio , pressure_ratio);    

    consv[3] = energy*rho + 0.5*rho*(prim[1]*prim[1] + prim[2]*prim[2]);
    return consv;
}
array ConservativeToPrimative(array consv){
    array prim;

    //find the simpler variables 
    prim[0] = consv[0]; //rho
    prim[1] = consv[1] / consv[0]; //u_x
    prim[2] = consv[2] / consv[0]; //u_y

    //find our index for density with binary search
    int point1 =0;
    int point2 =499;
    int half;
    int density_index;
    double rho = prim[0];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(rho > densities[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up the linear interpolation
    double density_top = densities[point2];
    int density_top_i = point2;
    double density_bottom = densities[point1];
    int density_bottom_i = point1;
    double density_ratio = (rho - density_bottom) / (density_top - density_bottom);
    if(density_bottom == density_top){density_ratio = 1;}

    //find our indices for pressure with binary search on the density-based linear interpolation of energies
    point1 =0;
    point2 =499;
    half;
    int pressure_index;
    double e = (consv[3] - 0.5*rho*(prim[1]*prim[1] + prim[2]*prim[2])) / rho; //internal energy

    //find the indices either side of our energy value in density top and bottom, they should be the same
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(e > (energies[half][density_bottom_i] + density_ratio*(energies[half][density_top_i] - energies[half][density_bottom_i]))){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }
    int pressure_top_i = point2;
    int pressure_bottom_i = point1;

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up our quadrant to search in
    double energies_NE = energies[pressure_top_i][density_top_i];
    double energies_SE = energies[pressure_bottom_i][density_top_i];
    double energies_SW = energies[pressure_bottom_i][density_bottom_i];
    double energies_NW = energies[pressure_top_i][density_bottom_i];

    //reverse bilinear interpolation
    double pressure_ratio = ReverseBilinearInterpolation(energies_NE , energies_SE , energies_SW , energies_NW , density_ratio , e); 
    
    //calculate the value of pressure
    double p = pressure_ratio * (pressures[pressure_top_i] - pressures[pressure_bottom_i]) + pressures[pressure_bottom_i];

    prim[3] = p;

    return prim;
}

//function to interpolate generally
double interpolate(array u, data_table dataTable){
    array prim = ConservativeToPrimative(u);

    //find our index for density with binary search
    int point1 =0;
    int point2 =499;
    int half;
    int density_index;
    double rho = u[0]; 
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(rho > densities[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}
    
    //set up the linear interpolation
    double density_top = densities[point2];
    int density_top_i = point2;
    double density_bottom = densities[point1];
    int density_bottom_i = point1;
    double density_ratio = (rho - density_bottom) / (density_top - density_bottom);
    if(density_bottom == density_top){density_ratio = 1;}
    

    //find our index for pressure also with linear bisection
    point1 =0;
    point2 =500;
    half;
    int pressure_index;
    double p = prim[2];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(p > pressures[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up linear interpolation
    double pressure_top = pressures[point2];
    int pressure_top_i = point2;
    double pressure_bottom = pressures[point1];
    int pressure_bottom_i = point1;
    double pressure_ratio = (p - pressure_bottom) / (pressure_top - pressure_bottom);
    if(pressure_bottom == pressure_top){pressure_ratio =1;}

    //then apply to the energies table
    double NE = dataTable[pressure_top_i][density_top_i];
    double SE = dataTable[pressure_bottom_i][density_top_i];
    double SW = dataTable[pressure_bottom_i][density_bottom_i];
    double NW = dataTable[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double value = BilinearInterpolation(NE , SE , SW , NW , density_ratio , pressure_ratio);  
    
    return value;
}

//function to find the sound speed with bilinear interpolation, assumming you know the conservative form
double SoundSpeed(array u){
    array prim = ConservativeToPrimative(u);

    //find our index for density with binary search
    int point1 =0;
    int point2 =499;
    int half;
    int density_index;
    double rho = u[0]; 
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(rho > densities[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}
    
    //set up the linear interpolation
    double density_top = densities[point2];
    int density_top_i = point2;
    double density_bottom = densities[point1];
    int density_bottom_i = point1;
    double density_ratio = (rho - density_bottom) / (density_top - density_bottom);
    if(density_bottom == density_top){density_ratio = 1;}
    

    //find our index for pressure also with linear bisection
    point1 =0;
    point2 =500;
    half;
    int pressure_index;
    double p = prim[3];
    while(std::abs(point1-point2)>1){
        half = floor((point1+point2)/2);
        if(p > pressures[half]){
            point1 = half;
        }
        else{
            point2 = half;
        }
    }

    //handle boundary terms
    if(point1 >= 499){point1 = 499;}
    if(point2 >= 499){point2 = 499;}
    if(point1 <=0){point1 = 0;}
    if(point2 <=0){point2 = 0;}

    //set up linear interpolation
    double pressure_top = pressures[point2];
    int pressure_top_i = point2;
    double pressure_bottom = pressures[point1];
    int pressure_bottom_i = point1;
    double pressure_ratio = (p - pressure_bottom) / (pressure_top - pressure_bottom);
    if(pressure_bottom == pressure_top){pressure_ratio =1;}

    //then apply to the energies table
    double ss_NE = sound_speeds[pressure_top_i][density_top_i];
    double ss_SE = sound_speeds[pressure_bottom_i][density_top_i];
    double ss_SW = sound_speeds[pressure_bottom_i][density_bottom_i];
    double ss_NW = sound_speeds[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double sound_speed = BilinearInterpolation(ss_NE , ss_SE , ss_SW , ss_NW , density_ratio , pressure_ratio);  
    
    return sound_speed;
}

//define the flux function but now we split into two functions fluxX and fluxY
array fluxX_def(std::array <double ,4> x){
    array fluxX;
    array prim = ConservativeToPrimative(x);
    double rho = prim[0];
    double v_x = prim[1];
    double v_y = prim[2];
    double p = prim[3];

    fluxX[0] = x[1];
    fluxX[1] = rho*v_x*v_x + p;
    fluxX[2] = rho * v_x * v_y;
    fluxX[3] = (x[3] + p)*v_x;

    return fluxX;
} //flux for the x-direction
array fluxY_def(std::array <double ,4> x){
    array fluxY;
    array prim = ConservativeToPrimative(x);
    double rho = prim[0];
    double v_x = prim[1];
    double v_y = prim[2];
    double p = prim[3];

    fluxY[0] = x[2];
    fluxY[1] = rho*v_x*v_y;
    fluxY[2] = rho * v_y * v_y + p;
    fluxY[3] = (x[3] + p)*v_y;

    return fluxY;
} //flux for the y-direction

array getFluxX(array x , array y , double dx , double dt){
    //find the flux at u_{i} and u_{i+1}

    array f_1 = fluxX_def(x); //f(ui)
    array f_2 = fluxX_def(y); //f(ui+1)

    //find u_{i+1/2}

    array uPlusHalf;
    for(int i=0; i<=3; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dx)*(f_2[i] - f_1[i]);
    }

    array RI_flux = fluxX_def(uPlusHalf); //richtmyer flux
    array LF_flux; 
    array FORCE_flux;
    
    for(int i=0; i<=3; ++i){
        LF_flux[i] = (0.5*(dx/dt))*(x[i] - y[i]) + 0.5*(f_1[i] + f_2[i]); // lax friedrichs flux
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]); //force flux
    }

    return FORCE_flux;

    
} //FORCE in x
array getFluxY(array x , array y , double dy , double dt){
    //find the flux at u_{i} and u_{i+1}

    array f_1 = fluxY_def(x); //f(ui)
    array f_2 = fluxY_def(y); //f(ui+1)

    //find u_{i+1/2}

    array uPlusHalf;
    for(int i=0; i<=3; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dy)*(f_2[i] - f_1[i]);
    }

    array RI_flux = fluxY_def(uPlusHalf); //richtmyer flux
    array LF_flux; 
    array FORCE_flux;
    
    for(int i=0; i<=3; ++i){
        LF_flux[i] = (0.5*dy/dt)*(x[i] - y[i]) + 0.5*(f_1[i] + f_2[i]); // lax friedrichs flux
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]); //force flux
    }

    return FORCE_flux;
    
} //FORCE in y

double ComputeTimeStep(const std::vector<std::vector<std::array <double , 4>>>& u , double C , double dx , double dy){
    double dspace = dx;
    if(dy < dx){
        dspace = dy;
    }

    double maxSpeed = 0.0;


    for(int i=2 ; i<=u.size()-2 ; ++i){
        for(int j=2; j<=u[0].size()-2; ++j){
            array prim = ConservativeToPrimative(u[i][j]);
            double rho = prim[0];
            double v_x = prim[1];
            double v_y = prim[2];
            double pressure = prim[3];
            double speed = std::sqrt(v_x*v_x + v_y*v_y);

            double sound_speed = SoundSpeed(u[i][j]);
            double new_speed = speed + sound_speed;

            if(new_speed > maxSpeed){
                maxSpeed = new_speed;
            }
            
        }
    }  
    
    double dt = C * dspace / maxSpeed;
    return dt;

} 

// Function to apply boundary conditions
void applyBoundaryConditions(std::vector<std::vector<std::array<double, 4>>>& u, int nxCells, int nyCells) {
    // ------ TRANSMISSIVE --------

    // // Left and Right boundaries
    // for (int j = 0; j < nyCells + 4; ++j) {
    //     u[0][j] = u[2][j];      // Left boundary: copy from first interior cell
    //     u[1][j] = u[2][j];      // Left ghost cell
    //     u[nxCells + 2][j] = u[nxCells + 1][j];  // Right ghost cell
    //     u[nxCells + 3][j] = u[nxCells + 1][j];  // Right boundary
    // }

    // // Bottom and Top boundaries
    // for (int i = 0; i < nxCells + 4; ++i) {
    //     u[i][0] = u[i][2];      // Bottom boundary: copy from first interior cell
    //     u[i][1] = u[i][2];      // Bottom ghost cell
    //     u[i][nyCells + 2] = u[i][nyCells + 1];  // Top ghost cell
    //     u[i][nyCells + 3] = u[i][nyCells + 1];  // Top boundary
    // }


    // -------- PERIODIC --------

    //left and right
    // for (int j = 0; j < nyCells + 4; ++j) {
    //     u[0][j] = u[nxCells][j]; 
    //     u[1][j] = u[nxCells+1][j];      // Left boundary
    //     u[nxCells + 2][j] = u[2][j]; 
    //     u[nxCells + 3][j] = u[3][j];  // Right boundary
    // }

    //bottom and top
    // for (int i = 0; i < nxCells + 4; ++i) {
    //     u[i][0] = u[i][nyCells];      // Bottom boundary
    //     u[i][1] = u[i][nyCells+1];    
    //     u[i][nyCells + 2] = u[i][2];  // Top boundary
    //     u[i][nyCells + 3] = u[i][3];  
    // }


    // ------ REFLECTIVE -------

    //left and right
    for (int j = 0; j < nyCells + 4; ++j) {//reflect in u_y 
        u[0][j] = u[2][j];
        u[1][j] = u[3][j];        // Bottom boundary
        u[nyCells + 2][j] = u[nyCells][j];
        u[nyCells + 3][j] = u[nyCells + 1][j];  // Top boundary
        u[0][j][2] = -u[3][j][2]; 
        u[1][j][2] = -u[2][j][2]; 
        u[nyCells + 2][j][2] = -u[nyCells + 1][j][2]; 
        u[nyCells + 3][j][2] = -u[nyCells][j][2]; 
    }

    // // Bottom and Top boundaries
    for (int i = 0; i < nxCells + 4; ++i) {//reflect in u_y
        u[i][0] = u[i][2];
        u[i][1] = u[i][3];        // Bottom boundary
        u[i][nyCells + 2] = u[i][nyCells];
        u[i][nyCells + 3] = u[i][nyCells + 1];  // Top boundary
        u[i][0][2] = -u[i][3][2]; 
        u[i][1][2] = -u[i][2][2]; 
        u[i][nyCells + 2][2] = -u[i][nyCells + 1][2]; 
        u[i][nyCells + 3][2] = -u[i][nyCells][2]; 
    }
}

// Update source terms
big_array SourceTermUpdate(big_array u , double x0,double dx, double y0, double dy, double dt){
    big_array update;
    update.resize(u.size(), std::vector<std::array<double, 4> >(u[0].size()));
    double alpha = 1.0;
    for (int j=0; j < u[0].size(); j++){
        for(int i = 0; i < u.size(); i++) { 
        
            double r = x0 + (i-0.5) * dx;
            double z = y0 + (j-0.5) * dy;
            double rho = u[i][j][0];
            double mom_r = u[i][j][1];
            double mom_z = u[i][j][2];
            double E = u[i][j][3];

            array prim = ConservativeToPrimative(u[i][j]);
            double v = prim[1];
            double w = prim[2];
            double p = prim[3];

            // Source terms
            double S_rho = -alpha * rho * v / r;
            double S_mom_r = -alpha * rho * v * v / r;
            double S_mom_z = -alpha * rho * v * w / r;
            double S_E = -alpha * (E + p) * v / r;

            // Update using RK2 for source term only:

            // Stage 1
            array u_stage1;
            u_stage1[0] = rho + dt * S_rho;
            u_stage1[1] = mom_r + dt * S_mom_r;
            u_stage1[2] = mom_z + dt * S_mom_z;
            u_stage1[3] = E + dt * S_E;

            // Recompute primitives for stage 2
            array prim_stage1 = ConservativeToPrimative(u_stage1);
            double v1 = prim_stage1[1];
            double w1 = prim_stage1[2];
            double p1 = prim_stage1[3];

            // Stage 2 source terms
            double S_rho_2 = -alpha * u_stage1[0] * v1 / r;
            double S_mom_r_2 = -alpha * u_stage1[0] * v1 * v1 / r;
            double S_mom_z_2 = -alpha * u_stage1[0] * w1 * v1 / r;
            double S_E_2 = -alpha * (u_stage1[3] + p1) * v1 / r;

            // Final update
            update[i][j][0] = rho + 0.5 * dt * (S_rho + S_rho_2);
            update[i][j][1] = mom_r + 0.5 * dt * (S_mom_r + S_mom_r_2);
            update[i][j][2] = mom_z + 0.5 * dt * (S_mom_z + S_mom_z_2);
            update[i][j][3] = E + 0.5 * dt * (S_E + S_E_2);
            
        }
    }
    return update;
}

std::vector<std::vector<std::array<double, 4> > > XthenY(std::vector<std::vector<std::array<double, 4> > > u , double dx , double dy , double dt , int nxCells , int nyCells , double x0 ,double x1 , double y0 ,double y1 ,double tStart , double tStop , double C,double omega ){
    applyBoundaryConditions(u, nxCells , nyCells);
    std::vector<std::vector<std::array<double, 4> > > uPlus1;
    uPlus1.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up uPlus1

    std::vector<std::vector<std::array<double, 4> > > fluxX;
    fluxX.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxX

    std::vector<std::vector<std::array<double, 4> > > fluxY;
    fluxY.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxY

    std::vector<std::vector<std::array<double, 4> > > uPlus1Bar;
    uPlus1Bar.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //uPlus1Bar

    std::vector<std::vector<std::array<double, 4> > > uBarL;
    uBarL.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarR;
    uBarR.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarHalfL;
    uBarHalfL.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarHalfR;
    uBarHalfR.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    

        //calculate the flux in the x-direction
        

        //start with uBar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    double DeltaPlus = u[i+1][j][k] - u[i][j][k];
                    double DeltaMinus = u[i][j][k] - u[i-1][j][k];
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    uBarL[i][j][k] = u[i][j][k] - 0.5 * xi * Delta;
                    uBarR[i][j][k] = u[i][j][k] + 0.5 * xi * Delta;
                    
                }

            }
        }
        //then uBarHalf
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j])[k]-fluxX_def(uBarL[i][j])[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j])[k]-fluxX_def(uBarL[i][j])[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in x direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxX[i][j] = getFluxX(uBarHalfR[i][j], uBarHalfL[i+1][j], dx, dt);
            }
        }
        //from this we get our intermediate ubar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1Bar[i][j][v] = u[i][j][v] - (dt/dx)*(fluxX[i][j][v]-fluxX[i-1][j][v]);
                }
            }
        }
        // apply boundary conditions to uPlus1Bar
        applyBoundaryConditions(uPlus1Bar , nxCells , nyCells);


        //now move in the y direction


        //start with uBar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    double DeltaPlus = uPlus1Bar[i][j+1][k] - uPlus1Bar[i][j][k];
                    double DeltaMinus = uPlus1Bar[i][j][k] - uPlus1Bar[i][j-1][k];
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    uBarL[i][j][k] = uPlus1Bar[i][j][k] - 0.5 * xi * Delta;  
                    uBarR[i][j][k] = uPlus1Bar[i][j][k] + 0.5 * xi * Delta;  
                    
                }

            }
        }
        //then uBarHalf
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j])[k]-fluxY_def(uBarL[i][j])[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j])[k]-fluxY_def(uBarL[i][j])[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxY(uBarHalfR[i][j], uBarHalfL[i][j+1], dy, dt);
            }
        }


        //update the results!
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1[i][j][v] = uPlus1Bar[i][j][v] - (dt/dy)*(fluxY[i][j][v]-fluxY[i][j-1][v]);
                }
            }
        }
    
    return uPlus1;
}

std::vector<std::vector<std::array<double, 4> > > YthenX(std::vector<std::vector<std::array<double, 4> > > u , double dx , double dy , double dt , int nxCells , int nyCells , double x0 ,double x1 , double y0 ,double y1 ,double tStart , double tStop , double C,double omega ){
    applyBoundaryConditions(u, nxCells , nyCells);
    std::vector<std::vector<std::array<double, 4> > > uPlus1;
    uPlus1.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up uPlus1

    std::vector<std::vector<std::array<double, 4> > > fluxX;
    fluxX.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxX

    std::vector<std::vector<std::array<double, 4> > > fluxY;
    fluxY.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //fluxY

    std::vector<std::vector<std::array<double, 4> > > uPlus1Bar;
    uPlus1Bar.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //uPlus1Bar

    std::vector<std::vector<std::array<double, 4> > > uBarL;
    uBarL.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarR;
    uBarR.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarHalfL;
    uBarHalfL.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    std::vector<std::vector<std::array<double, 4> > > uBarHalfR;
    uBarHalfR.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4));

    

        //calculate the flux in the x-direction
        

        //start with uBar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    double DeltaPlus = u[i][j+1][k] - u[i][j][k];
                    double DeltaMinus = u[i][j][k] - u[i][j-1][k];
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    uBarL[i][j][k] = u[i][j][k] - 0.5 * xi * Delta;
                    uBarR[i][j][k] = u[i][j][k] + 0.5 * xi * Delta;
                    
                }

            }
        }
        //then uBarHalf
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j] )[k]-fluxY_def(uBarL[i][j] )[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dy)*(fluxY_def(uBarR[i][j] )[k]-fluxY_def(uBarL[i][j] )[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxY(uBarHalfR[i][j], uBarHalfL[i][j+1], dy, dt);
            }
        }
        //from this we get our intermediate ubar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1Bar[i][j][v] = u[i][j][v] - (dt/dy)*(fluxY[i][j][v]-fluxY[i][j-1][v]);
                }
            }
        }
        // apply boundary conditions to uPlus1Bar
        applyBoundaryConditions(uPlus1Bar , nxCells , nyCells);


        //now move in the x direction


        //start with uBar
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    double DeltaPlus = uPlus1Bar[i+1][j][k] - uPlus1Bar[i][j][k];
                    double DeltaMinus = uPlus1Bar[i][j][k] - uPlus1Bar[i-1][j][k];
                    double r = DeltaMinus / (DeltaPlus  + 1e-8);
                    
                    double xi_L = 2.0*r/(1+r);
                    double xi_R = 2.0/(1+r);
                    double xi;
                    double Delta = 0.5*(1+omega)*DeltaMinus + 0.5*(1-omega)*DeltaPlus;
                    
                    if(r<=0){ xi=0;}
                    else if(r>0 && r<=1){ xi=r;}
                    else{ 
                        xi=std::fmin(1, xi_R);
                        
                    }

                    uBarL[i][j][k] = uPlus1Bar[i][j][k] - 0.5 * xi * Delta;  
                    uBarR[i][j][k] = uPlus1Bar[i][j][k] + 0.5 * xi * Delta;  
                    
                }

            }
        }
        //then uBarHalf
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    uBarHalfL[i][j][k] = uBarL[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j] )[k]-fluxX_def(uBarL[i][j] )[k]);
                    uBarHalfR[i][j][k] = uBarR[i][j][k] - 0.5*(dt/dx)*(fluxX_def(uBarR[i][j] )[k]-fluxX_def(uBarL[i][j] )[k]);
                }
            }
        }
        applyBoundaryConditions(uBarHalfL, nxCells , nyCells);
        applyBoundaryConditions(uBarHalfR, nxCells , nyCells);
        //now we have flux in y direction
        for (int j = 1; j < nyCells + 3; j++) {
            for (int i = 1; i < nxCells + 3; i++) { 
                fluxY[i][j] = getFluxX(uBarHalfR[i][j], uBarHalfL[i+1][j], dx, dt);
            }
        }


        //update the results!
        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int v = 0; v <= 3; v++) {
                    uPlus1[i][j][v] = uPlus1Bar[i][j][v] - (dt/dx)*(fluxX[i][j][v]-fluxX[i-1][j][v]);
                }
            }
        }
    
    return uPlus1;
}

void SaveData(const big_array& results, double x0, double dx, double y0, double dy, int nxCells, int nyCells, int index, int nTheta = 200) {
    std::string filename = "SLIC" + std::to_string(index) + ".dat";
    std::ofstream output(filename);
    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            double x = x0 + (i - 1)*dx;
            double y = y0 + (j - 1)*dy;
            output << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] << std::endl;
        }
        output<<std::endl;
    }
    output.close();
}
//flatten the (i,j) into the index of a vector of size N (there will be Nz of each i with each j added on)
int vectoridx(int i, int j, int Nz){
    return i * Nz + j;
}

//solve the sparse matrix linear system for laplace's equation
Eigen::VectorXd laplaceSolver(big_array& u, double t, double x0, double dx, double y0, double dy, double Nr, double Nz){
    double dr = dx;
    double dz = dy;
    double dr2 = dr * dr;
    double dz2 = dz * dz;
    int N = Nr * Nz;
    double r0 = x0;
    double z0 = y0;
    int electrode = std::floor(Nr / 3); // set the electrode to be at 1/3 of the domain

    Eigen::SparseMatrix<double> M(N,N);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
    std::vector<Eigen::Triplet<double>> coefficients;

    //set up matrix
    for(int i =0; i < Nr; ++i){
        double r = r0 + (i - 1.5)* dr;

        for(int j =0; j < Nz; ++j){
            int k = vectoridx(i,j, Nz);
            double sigma = interpolate(u[i][j], electrical_conductivity);

            if(i < electrode && j == Nz - 1){
                coefficients.emplace_back(k,k,1.0);
                b[k] = -(1/sigma) * (I0 * (std::exp(-alpha * t) * std::sin))/(PI*R0*R0);
                continue;
            }

            if(i > electrode && j == Nz - 1){
                // Neumann at the rest of the top
                coefficients.emplace_back(k, vectoridx(i, j, Nz), 1.0);
                coefficients.emplace_back(k, vectoridx(i, j-1, Nz), -1.0);
                //coefficients.emplace_back(k,k,1.0);
                continue;
            }

            if(i == 0){
                // Neumann at r = r0 boundary:
                coefficients.emplace_back(k, vectoridx(i, j, Nz), -1.0);
                coefficients.emplace_back(k, vectoridx(i+1, j, Nz), 1.0);
                //coefficients.emplace_back(k,k,1.0);
                continue;
            }

            if(i == Nr - 1){
                // Neumann at outer r boundary
                coefficients.emplace_back(k, vectoridx(i, j, Nz), 1.0);
                coefficients.emplace_back(k, vectoridx(i-1, j, Nz), -1.0);
                //coefficients.emplace_back(k,k,1.0);
                continue;
            }

            if(j == 0){
                // Neumann at z = z0 boundary
                coefficients.emplace_back(k, vectoridx(i, j, Nz), -1.0);
                coefficients.emplace_back(k, vectoridx(i, j+1, Nz), 1.0);
                //coefficients.emplace_back(k,k,1.0);
                continue;
            }
            

            //coefficients for finite difference scheme for poisson equation
            double rPlus = r + dr / 2.0;
            double rMinus = r - dr / 2.0;
            double centre = -2.0 / dr2 - 2.0 / dz2;
            double coeffPlusr = (1.0 + dr / (2.0*r)) / dr2;
            double coeffMinusr = (1.0 - dr / (2.0*r)) / dr2;
            double coeffz = 1.0 / dz2;

            // Fill matrix entries
            coefficients.emplace_back(k, vectoridx(i, j, Nz), centre);
            coefficients.emplace_back(k, vectoridx(i + 1, j, Nz), coeffPlusr);
            coefficients.emplace_back(k, vectoridx(i - 1, j, Nz), coeffMinusr);
            coefficients.emplace_back(k, vectoridx(i, j + 1, Nz), coeffz);
            coefficients.emplace_back(k, vectoridx(i, j - 1, Nz), coeffz);

            
        }
    }
    //build a sparse matrix
    M.setFromTriplets(coefficients.begin(), coefficients.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::VectorXd phi = solver.solve(b);

    return phi;
}

//set up current distribution
std::array<double,2> current(double r, double z, const Eigen::VectorXd& phi, double t, double x0, double dx, double y0, double dy, double Nr, double Nz, double gamma = 10847100){
    int i = static_cast<int>(std::round((r - x0)/dx + 1.5));
    int j = static_cast<int>(std::round((z - y0)/dy + 1.5));

    // Ensure i and j are within valid bounds (1..Nr-2 and 1..Nz-2 for central differences)
    if (i <= 0){i = 1;}
    if (i >= Nr - 1){i = Nr - 2;}
    if (j <= 0){ j = 1;}
    if (j >= Nz - 1){ j = Nz - 2;}

    int k_ip = vectoridx(i+1, j, Nz);
    int k_im = vectoridx(i-1, j, Nz);
    int k_jp = vectoridx(i, j+1, Nz);
    int k_jm = vectoridx(i, j-1, Nz);

    double Jr = - 312*(phi[k_ip] - phi[k_im]) / (2*dx);
    double Jz = - 312*(phi[k_jp] - phi[k_jm]) / (2*dy);



    // double I = I0 * (std::exp(-alpha * t) - std::exp(-beta * t)) * ( 1 - std::exp(-gamma * t)) * ( 1 - std::exp(-gamma * t));
    std::array<double,2> J;
    // J[1] = - (I/PI*R0*R0) * std::exp(-(r/R0)*(r/R0));
    // J[0] = 0;

    J[0] = Jr;
    J[1] = Jz;
    return J;
}

//solve the sparse matrix linear system
std::tuple<Eigen::VectorXd, Eigen::VectorXd> poissonSolverR(double t, double x0, double dx, double y0, double dy, double Nr, double Nz){
    double dr = dx;
    double dz = dy;
    double dr2 = dr * dr;
    double dz2 = dz * dz;
    int N = Nr * Nz;
    double r0 = x0;
    double z0 = y0;
    Eigen::VectorXd phi = laplaceSolver(t,x0,dx,y0,dy,Nr,Nz);

    Eigen::SparseMatrix<double> M(N,N);
    Eigen::VectorXd bR = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd bZ = Eigen::VectorXd::Zero(N);
    std::vector<Eigen::Triplet<double>> coefficients;

    //set up matrix
    for(int i =0; i < Nr; ++i){
        double r = r0 + (i - 1.5)* dr;

        for(int j =0; j < Nz; ++j){
            int k = vectoridx(i,j, Nz);
            double z = z0 + (j - 1.5) * dz;

            //RHS
            bR[k] = -current(r, z, phi, t, x0, dx, y0, dy, Nr, Nz)[0];
            bZ[k] = -current(r, z, phi, t, x0, dx, y0, dy, Nr, Nz)[1];

            //boundary conditions
            if(j == Nz - 1){
                //dirichlet conditions
                coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }
            if(i == 0){
                // Neumann at r = r0 boundary:
                coefficients.emplace_back(k, vectoridx(i, j, Nz), -1.0);
                coefficients.emplace_back(k, vectoridx(i+1, j, Nz), 1.0);
                 //coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }

            if(i == Nr - 1){
                // Neumann at outer r boundary
                coefficients.emplace_back(k, vectoridx(i, j, Nz), 1.0);
                coefficients.emplace_back(k, vectoridx(i-1, j, Nz), -1.0);
                 //coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }

            if(j == 0){
                //dirichlet conditions
                coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }
        
            //coefficients for finite difference scheme for poisson equation
            double rPlus = r + dr / 2.0;
            double rMinus = r - dr / 2.0;
            double centre = -2.0 / dr2 - 2.0 / dz2;
            double coeffPlusr = (1.0 + dr / (2.0*r)) / dr2;
            double coeffMinusr = (1.0 - dr / (2.0*r)) / dr2;
            double coeffz = 1.0 / dz2;

            // Fill matrix entries
            coefficients.emplace_back(k, vectoridx(i, j, Nz), centre);
            coefficients.emplace_back(k, vectoridx(i + 1, j, Nz), coeffPlusr);
            coefficients.emplace_back(k, vectoridx(i - 1, j, Nz), coeffMinusr);
            coefficients.emplace_back(k, vectoridx(i, j + 1, Nz), coeffz);
            coefficients.emplace_back(k, vectoridx(i, j - 1, Nz), coeffz);
        
            
        }
    }

    //build a sparse matrix
    M.setFromTriplets(coefficients.begin(), coefficients.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::VectorXd AR = solver.solve(bR);

    return {AR,bR};
}

//solve the sparse matrix linear system
std::tuple<Eigen::VectorXd, Eigen::VectorXd> poissonSolverZ(double t, double x0, double dx, double y0, double dy, double Nr, double Nz){
    double dr = dx;
    double dz = dy;
    double dr2 = dr * dr;
    double dz2 = dz * dz;
    int N = Nr * Nz;
    double r0 = x0;
    double z0 = y0;
    Eigen::VectorXd phi = laplaceSolver(t,x0,dx,y0,dy,Nr,Nz);

    Eigen::SparseMatrix<double> M(N,N);
    Eigen::VectorXd bR = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd bZ = Eigen::VectorXd::Zero(N);
    std::vector<Eigen::Triplet<double>> coefficients;

    //set up matrix
    for(int i =0; i < Nr; ++i){
        double r = r0 + (i - 1.5)* dr;

        for(int j =0; j < Nz; ++j){
            int k = vectoridx(i,j, Nz);
            double z = z0 + (j - 1.5) * dz;

            //RHS
            bR[k] = -current(r, z, phi, t, x0, dx, y0, dy, Nr, Nz)[0];
            bZ[k] = -current(r, z, phi, t, x0, dx, y0, dy, Nr, Nz)[1];

            //boundary conditions
            if(j == Nz - 1){
                // Neumann 
                coefficients.emplace_back(k, vectoridx(i, j, Nz), 1.0);
                coefficients.emplace_back(k, vectoridx(i, j-1, Nz), -1.0);
                 //coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }
            if(i == 0){
                // Neumann 
                coefficients.emplace_back(k, vectoridx(i, j, Nz), 1.0);
                coefficients.emplace_back(k, vectoridx(i+1, j, Nz), -1.0);
                //coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }

            if(i == Nr - 1){
                //dirichlet conditions
                coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }

            if(j == 0){
                // Neumann 
                coefficients.emplace_back(k, vectoridx(i, j, Nz), 1.0);
                coefficients.emplace_back(k, vectoridx(i, j+1, Nz), -1.0);
                 //coefficients.emplace_back(k,k,1.0);
                bR[k] = 0.0;
                bZ[k] = 0.0;
                continue;
            }
          
            //coefficients for finite difference scheme for poisson equation
            double rPlus = r + dr / 2.0;
            double rMinus = r - dr / 2.0;
            double centre = -2.0 / dr2 - 2.0 / dz2;
            double coeffPlusr = (1.0 + dr / (2.0*r)) / dr2;
            double coeffMinusr = (1.0 - dr / (2.0*r)) / dr2;
            double coeffz = 1.0 / dz2;

            // Fill matrix entries
            coefficients.emplace_back(k, vectoridx(i, j, Nz), centre);
            coefficients.emplace_back(k, vectoridx(i + 1, j, Nz), coeffPlusr);
            coefficients.emplace_back(k, vectoridx(i - 1, j, Nz), coeffMinusr);
            coefficients.emplace_back(k, vectoridx(i, j + 1, Nz), coeffz);
            coefficients.emplace_back(k, vectoridx(i, j - 1, Nz), coeffz);
        
        }
    }

    //build a sparse matrix
    M.setFromTriplets(coefficients.begin(), coefficients.end());

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::VectorXd AZ = solver.solve(bZ);

    return {AZ,bZ};
}

//function to now retrive the magnetic feild from the magnetic potential
std::vector<double> magneticfield1(double t, double x0, double dx, double y0, double dy, int nxCells, int nyCells){
    auto [Ar,Jr] = poissonSolverR(t, x0, dx,y0, dy,nxCells, nyCells);
    auto [Az,Jz] = poissonSolverZ(t, x0, dx,y0, dy,nxCells, nyCells);

    //need to find first derviatives in both for the magnetic field
    std::vector<double> Br(Ar.size());
    //central differences
    for (int i=1; i<Ar.size()-1; i++){
        Br[i] = - (Ar[i+1] - Ar[i-1]) / (2.0 * dy);
    }

    //endpoints
    Br[0] = -(Ar[1] - Ar[0]) / dy;
    Br[Br.size() - 1] = -(Ar[Ar.size() - 1] - Ar[Ar.size() - 2]) / dy;

    std::vector<double> Bz(Az.size());
    //central differences
    for (int i=1; i<Az.size()-1; i++){
        Bz[i] = - (Az[i+1] - Az[i-1]) / (2.0 * dx);
    }

    //endpoints
    Bz[0] = -(Az[1] - Az[0]) / dx;
    Bz[Bz.size() - 1] = -(Az[Az.size() - 1] - Az[Az.size() - 2]) / dx;

    //final magnetic field
    std::vector<double> B(Ar.size());
    for ( int i =0; i<Ar.size(); i++){
        B[i] = Br[i] - Bz[i];
    }

    return B;
}

std::vector<double> magneticfield(double t, double x0, double dx, double y0, double dy, int nxCells, int nyCells){
    auto [Ar,Jr] = poissonSolverR(t, x0, dx,y0, dy,nxCells, nyCells);
    auto [Az,Jz] = poissonSolverZ(t, x0, dx,y0, dy,nxCells, nyCells);

    std::vector<double> Btheta(nxCells * nyCells, 0.0);

    // Interior points (central differences)
    for (int i = 1; i < nxCells - 1; i++) {
        for (int j = 1; j < nyCells - 1; j++) {
            int k = vectoridx(i, j, nyCells);

            double dAr_dz = (Ar[vectoridx(i, j+1, nyCells)] - Ar[vectoridx(i, j-1, nyCells)]) / (2.0 * dy);
            double dAz_dr = (Az[vectoridx(i+1, j, nyCells)] - Az[vectoridx(i-1, j, nyCells)]) / (2.0 * dx);

            Btheta[k] = dAr_dz - dAz_dr;
        }
    }

    // Boundaries - use one-sided finite differences

    // Left edge i = 0
    for (int j = 0; j < nyCells; j++) {
        int k = vectoridx(0, j, nyCells);
        int k_right = vectoridx(1, j, nyCells);
        int k_up = vectoridx(0, std::min(j+1, nyCells-1), nyCells);
        int k_down = vectoridx(0, std::max(j-1, 0), nyCells);

        double dAr_dz = (Ar[k_up] - Ar[k_down]) / (2.0 * dy);
        double dAz_dr = (Az[k_right] - Az[k]) / dx;  // forward difference

        Btheta[k] = dAr_dz - dAz_dr;
    }

    // Right edge i = nxCells - 1
    for (int j = 0; j < nyCells; j++) {
        int k = vectoridx(nxCells - 1, j, nyCells);
        int k_left = vectoridx(nxCells - 2, j, nyCells);
        int k_up = vectoridx(nxCells - 1, std::min(j+1, nyCells-1), nyCells);
        int k_down = vectoridx(nxCells - 1, std::max(j-1, 0), nyCells);

        double dAr_dz = (Ar[k_up] - Ar[k_down]) / (2.0 * dy);
        double dAz_dr = (Az[k] - Az[k_left]) / dx;  // backward difference

        Btheta[k] = dAr_dz - dAz_dr;
    }

    // Bottom edge j = 0
    for (int i = 0; i < nxCells; i++) {
        int k = vectoridx(i, 0, nyCells);
        int k_right = vectoridx(std::min(i+1, nxCells-1), 0, nyCells);
        int k_left = vectoridx(std::max(i-1, 0), 0, nyCells);
        int k_up = vectoridx(i, 1, nyCells);

        double dAr_dz = (Ar[k_up] - Ar[k]) / dy;  // forward difference
        double dAz_dr = (Az[k_right] - Az[k_left]) / (2.0 * dx);

        Btheta[k] = dAr_dz - dAz_dr;
    }

    // Top edge j = nyCells - 1
    for (int i = 0; i < nxCells; i++) {
        int k = vectoridx(i, nyCells - 1, nyCells);
        int k_right = vectoridx(std::min(i+1, nxCells-1), nyCells - 1, nyCells);
        int k_left = vectoridx(std::max(i-1, 0), nyCells - 1, nyCells);
        int k_down = vectoridx(i, nyCells - 2, nyCells);

        double dAr_dz = (Ar[k] - Ar[k_down]) / dy;  // backward difference
        double dAz_dr = (Az[k_right] - Az[k_left]) / (2.0 * dx);

        Btheta[k] = dAr_dz - dAz_dr;
    }
    return Btheta;
}


void momentumUpdate(big_array& u, double x0, double dx, double y0, double dy, double t, double dt, int nxCells, int nyCells, int max_iter = 10, double tol = 1e-6){
    auto [Az,Jz] = poissonSolverZ(t, x0, dx,y0, dy,nxCells, nyCells);
    auto [Ar,Jr] = poissonSolverR(t, x0, dx,y0, dy,nxCells, nyCells);

    std::vector<double> B(u.size());
    B = magneticfield(t,x0,dx,y0,dy, nxCells, nyCells);

    big_array update = u;
    std::vector<double> crossR(u.size());
    std::vector<double> crossZ(u.size());

    // set up subcycling
    int sub_steps = 10;
    double dt_sub = dt / sub_steps;

    for(int step = 0; step < sub_steps; ++step){

        for (int i =0; i<crossR.size(); i++){
            // J_z * B_{theta} (in -r direction)
            crossR[i] = -Jz[i] * B[i]; 
        }
        for (int i =0; i<crossZ.size(); i++){
            // J_r * B_{theta} (in z direction)
            crossZ[i] = Jr[i] * B[i]; 
        }

        for(int j =0; j<u[0].size(); j++){ // modelling the z direction
            for(int i = 0; i < u.size(); i++) { // modelling the r direction
            
                u[i][j][1] = u[i][j][1] + dt_sub * crossR[i];
                u[i][j][2] = u[i][j][2] + dt_sub * crossZ[i];
            }
        }

    }
}

void energyUpdate(big_array& u, double x0, double dx, double y0, double dy, double t, double dt, int nxCells, int nyCells, int max_iter = 10, double tol = 1e-6) {
    auto [Az,Jz] = poissonSolverZ(t, x0, dx,y0, dy,nxCells, nyCells);
    auto [Ar,Jr] = poissonSolverR(t, x0, dx,y0, dy,nxCells, nyCells);

    std::vector<double> B(u.size());
    B = magneticfield(t,x0,dx,y0,dy, nxCells, nyCells);

    big_array update = u;
    std::vector<double> crossR(u.size());
    std::vector<double> crossZ(u.size());

    // set up subcycling
    int sub_steps = 10;
    double dt_sub = dt / sub_steps;

    for(int step = 0; step < sub_steps; ++step){

        for (int i =0; i<crossR.size(); i++){
            // J_z * B_{theta} (in -r direction)
            crossR[i] = -Jz[i] * B[i]; 
        }
        for (int i =0; i<crossZ.size(); i++){
            // J_r * B_{theta} (in z direction)
            crossZ[i] = Jr[i] * B[i]; 
        }

        for(int j =0; j<u[0].size(); j++){ // modelling the z direction
            for(int i = 0; i < u.size(); i++) { // modelling the r direction
            
                double v_r = u[i][j][1] / u[i][j][0];
                double v_z = u[i][j][2] / u[i][j][0];

                double e1 = v_r * crossR[i] + v_z * crossZ[j];
                u[i][j][3] = u[i][j][3] + dt_sub * e1;
            }
        }

    }

}

void thermalSourceTerm_Newton(big_array& u, double dt, double t, double dx, double x0, double dy, double y0, int nxCells, int nyCells) {
    auto [Az,Jz] = poissonSolverZ(t, x0, dx,y0, dy,nxCells, nyCells);
    auto [Ar,Jr] = poissonSolverR(t, x0, dx,y0, dy,nxCells, nyCells);

    int max_iter = 20;
    double tol = 1e-16;
    // auto [A, J] = SolvePotential(t+dt, x0, dx, nCells);  // use t + dt for implicit

    double kappa = 60;
    double stefan_boltzmann = 5.67e-8;
    double epsilon = 1e-5;  // for finite difference

    for(int j = 0; j < u[0].size(); ++j){
        for (int i = 0; i < u.size(); ++i) {
            double E_old = u[i][j][3];  // E^n
            double E = E_old;        // initial guess E^0

            double rho = u[i][j][0];
            double momentum_r = u[i][j][1];
            double momentum_z = u[i][j][2];
            double velocity2 = (u[i][j][1] / u[i][j][0]) * (u[i][j][1] / u[i][j][0]) + (u[i][j][2] / u[i][j][0]) * (u[i][j][2] / u[i][j][0]);

            // Precompute J²
            double J2 = (Jz[i]*Jz[i]) + (Jr[i]*Jr[i]);

            // Compute constant mass-fraction-based heat of formation
            double heat_of_formation = 0.0;
            for (int j = 0; j < 19; ++j) {
                double mass_frac = interpolate(u[i][j], mass_fractions[j]);
                heat_of_formation += mass_frac * heats_of_formation[j];
            }

            int iter = 0;
            double diff = 1e9;

            while (iter < max_iter && diff > tol) {
                // T from energy guess
                double T = interpolate({rho, momentum_r, momentum_z, E}, temperatures);

                // σ from E guess
                double sigma = interpolate({rho, momentum_r, momentum_z, E}, electrical_conductivity);

                // Radiation loss term S_T
                double S_T = stefan_boltzmann * kappa * (std::pow(T, 4) - std::pow(T0, 4)) - rho * (heat_of_formation + 0.5 * velocity2);

                // f(E)
                double f = E - E_old - dt * ((1.0 / sigma) * J2 - S_T);

                // f'(E) via finite difference
                double E_eps = E + epsilon;

                // T_eps from E+ε
                double T_eps = interpolate({rho, momentum_r, momentum_z, E_eps}, temperatures);
                double sigma_eps = interpolate({rho, momentum_r, momentum_z, E_eps}, electrical_conductivity);
                double S_T_eps = stefan_boltzmann * kappa * (std::pow(T_eps, 4) - std::pow(T0, 4)) - rho * (heat_of_formation + 0.5 * velocity2);
                double f_eps = E_eps - E_old - dt * ((1.0 / sigma_eps) * J2 - S_T_eps);

                double df_dE = (f_eps - f) / epsilon;

                // Newton step
                double E_new = E - f / df_dE;
                diff = std::abs(E_new - E);
                E = E_new;
                ++iter;
            }
            // Update energy directly in the input array
            u[i][j][3] = E;
        }
    }
}


int main(){
    double tStart = 0.0;
    double tStop = 1.5e-4;
    double C = 0.8;
    double omega =0;
    int nxCells = 50;
    int nyCells = 50;
    double x0 = -1.0;
    double x1 = 1.0;
    double y0 = -1.0;
    double y1 = 1.0;

    std::vector<std::vector<std::array<double, 4> > > u;
    u.resize(nxCells+4, std::vector<std::array<double, 4> >(nyCells + 4)); //set up u

    double dx = (x1 - x0) / nxCells;
    double dy = (y1 - y0) / nyCells;

    //intial conditions
    for(int i = 0; i < u.size(); i++) { 
        for(int j = 0; j < u[0].size(); j++) {
            array prim;
            // coordinates
            double x = x0 + (i - 1.5) * dx;
            double y = y0 + (j - 1.5) * dy;

            
            if ( std::sqrt(x*x + y*y)<=0.4) {
                prim[0] = 1.0;
                prim[1] = 0.0*std::pow(10,2.5);
                prim[2] = 0.0*std::pow(10,2.5);
                prim[3] = 1.0*1e5;  // pressure (p)
            } 
            else{
                prim[0] = 0.125;
                prim[1] = 0.0*std::pow(10,2.5);
                prim[2] = 0.0*std::pow(10,2.5);
                prim[3] = 0.1*1e5;
            }

            // prim[0] = 1.225; // Density
            // prim[1] = 0*std::pow(10,2.5); // Velocity
            // prim[2] = 0;
            // prim[3] = 101325 + 2e6*std::exp(-((x*x+y*y)/R0)*((x*x+y*y)/R0));  // Pressure

            u[i][j] = PrimativeToConservative(prim);
        }
    }

    applyBoundaryConditions(u , nxCells , nyCells);
    double dt;
    double t = tStart;

    int counter =0;
    do{
        dt = ComputeTimeStep(u,C,dx,dy);
        t +=dt;

        // Update cylindrical source terms
        u = SourceTermUpdate(u,x0, dx,y0,dy,dt);

        std::cout << "t = "<< t<<" dt = "<< dt<< std::endl; 
        applyBoundaryConditions(u , nxCells , nyCells);

        std::vector<std::vector<std::array<double, 4> > > uPlus1x = XthenY(u ,  dx ,  dy ,  dt ,  nxCells ,  nyCells ,  x0 , x1 ,  y0 , y1 , tStart ,  tStop ,  C, omega);
        std::vector<std::vector<std::array<double, 4> > > uPlus1y = YthenX(u ,  dx ,  dy ,  dt ,  nxCells ,  nyCells ,  x0 , x1 ,  y0 , y1 , tStart ,  tStop ,  C, omega);

        for(int j = 2; j < nyCells+2; j++) {
            for(int i = 2; i < nxCells+2; i++) {
                for(int k=0 ; k<4 ; ++k){
                    u[i][j][k] = 0.5*(uPlus1x[i][j][k] + uPlus1y[i][j][k]);
                }
            }
        }
        applyBoundaryConditions(u , nxCells , nyCells);

        // Output data at specific time steps
        // while (t >= 1e-5 * counter) {
        //     std::vector<std::vector<std::array<double, 4> > > results;
        //     results.resize(nxCells+2, std::vector<std::array<double, 4> >(nyCells + 2)); //results
        
        //     for(int j = 1; j < nyCells+2; j++) { 
        //         for(int i = 1; i < nxCells+2; i++) {
        //             results[i][j] = ConservativeToPrimative(u[i][j]);
        //         }
        //     }

        //     SaveData(results, x0, dx, y0, dy, nxCells, nyCells, counter);
        //     std::cout << "Saved frame: " << counter << std::endl;
        //     counter += 1;
        // }

        //update thermal source terms
        thermalSourceTerm_Newton(u,dt,t,dx,x0,dy,y0,nxCells, nyCells);

        //update resistive source terms
        momentumUpdate(u, x0, dx,y0, dy, t, dt, nxCells, nyCells);
        energyUpdate(u, x0, dx,y0, dy, t, dt, nxCells, nyCells);
        
    }while(t<tStop);

    //now convert back to primitive
    std::vector<std::vector<std::array<double, 4> > > results;
    results.resize(nxCells+2, std::vector<std::array<double, 4> >(nyCells + 2)); //results
 
    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            results[i][j] = ConservativeToPrimative(u[i][j]);
        }
    }


    //output the results
    std::ofstream output("SLIC2d.dat");
    for(int j = 1; j < nyCells+2; j++) { 
        for(int i = 1; i < nxCells+2; i++) {
            double x = x0 + (i - 1)*dx;
            double y = y0 + (j - 1)*dy;
            output << x << " " << y << " " << results[i][j][0] << " " << results[i][j][1] << " " << results[i][j][2] << " " << results[i][j][3] << std::endl;
        }
        output<<std::endl;
    }

}