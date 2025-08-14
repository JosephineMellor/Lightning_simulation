#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <array>
#include <filesystem>
#include <tuple>

//defining some useful types
typedef std::array<double , 500> data_vec;
typedef std::vector<double> big_vector;
typedef std::vector<std::vector<double>> data_table;
typedef std::array<double,3> array;
typedef std::vector<std::array<double,3>> big_array;
typedef std::array<double , 19> species_vec;
typedef std::array<data_table, 19> species_tables;

const double PI = 3.141592653589793;
const double r0 = 2e-3; //2cm in meters
const double T0 = 298;
const double mu_0 = 4.0 * PI * 1e-7;

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
    double energies_NE = energies[pressure_top_i][density_top_i];
    double energies_SE = energies[pressure_bottom_i][density_top_i];
    double energies_SW = energies[pressure_bottom_i][density_bottom_i];
    double energies_NW = energies[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double energy = BilinearInterpolation(energies_NE , energies_SE , energies_SW , energies_NW , density_ratio , pressure_ratio);    

    consv[2] = energy*rho + 0.5*rho*(prim[1]*prim[1]);
    return consv;
}
array ConservativeToPrimative(array consv){
    array prim;

    //find the simpler variables 
    prim[0] = consv[0]; //rho
    prim[1] = consv[1] / consv[0]; //u_x

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
    double e = (consv[2] - 0.5*rho*(prim[1]*prim[1])) / rho; //internal energy

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

    prim[2] = p;

    return prim;
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
    double ss_NE = sound_speeds[pressure_top_i][density_top_i];
    double ss_SE = sound_speeds[pressure_bottom_i][density_top_i];
    double ss_SW = sound_speeds[pressure_bottom_i][density_bottom_i];
    double ss_NW = sound_speeds[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double sound_speed = BilinearInterpolation(ss_NE , ss_SE , ss_SW , ss_NW , density_ratio , pressure_ratio);  
    
    return sound_speed;
}

//function to find the temperature with bilinear interpolation, assumming you know the conservative form
double temperature(array u){
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
    double t_NE = temperatures[pressure_top_i][density_top_i];
    double t_SE = temperatures[pressure_bottom_i][density_top_i];
    double t_SW = temperatures[pressure_bottom_i][density_bottom_i];
    double t_NW = temperatures[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double temperature = BilinearInterpolation(t_NE , t_SE , t_SW , t_NW , density_ratio , pressure_ratio);  
    
    return temperature;
}

//function to find the thermal conductivity with bilinear interpolation, assumming you know the conservative form
double thermalConductivity(array u){
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
    double tc_NE = thermal_conductivity[pressure_top_i][density_top_i];
    double tc_SE = thermal_conductivity[pressure_bottom_i][density_top_i];
    double tc_SW = thermal_conductivity[pressure_bottom_i][density_bottom_i];
    double tc_NW = thermal_conductivity[pressure_top_i][density_bottom_i];

    //linear interpolate using our variables from before
    double thermalConductivity = BilinearInterpolation(tc_NE , tc_SE , tc_SW , tc_NW , density_ratio , pressure_ratio);  
    
    return thermalConductivity;
}

//function to save data as spherically symmetric
void SaveWrappedData(const std::vector<std::array<double, 3>>& results, double x0, double dx, int nCells, int index, int nTheta = 200) {
    std::string filename = "wrapped_" + std::to_string(index) + ".dat";
    std::ofstream out(filename);
    // std::ofstream out("wrapped.dat");
    for (int i = 0; i <= nCells+3; ++i) {
        double r = x0 + (i+0.5) * dx;
        double rho = results[i][0]; // Use density (or change to results[i][1] for momentum, etc.)

        for (int j = 0; j < nTheta; ++j) {
            double theta = 2.0 * PI * j / nTheta;
            double x = r * cos(theta);
            double y = r * sin(theta);
            out << x << " " << y << " " << rho << "\n";
        }
        out << "\n"; // Separate rings

        // out<< r << results[i][0]<<results[i][1]<<results[i][2]<<"\n";
    }
    out.close();
}




//function to save data as 1 dimensional
void SaveUnwrappedData(const std::vector<std::array<double, 3>>& u, double x0, double dx, int nCells, int index, int nTheta = 200) {
    big_array results(u.size());
    for(int i =0; i<u.size(); i++){
        results[i] = ConservativeToPrimative(u[i]);
    }
    
    std::string filename = "wrapped_" + std::to_string(index) + ".dat";
    std::ofstream out(filename);
    // std::ofstream out("wrapped.dat");
    for (int i = 0; i <= nCells+3; ++i) {
        double r = x0 + (i+0.5) * dx;
        double rho = results[i][0]; // Use density (or change to results[i][1] for momentum, etc.)

        double x = x0 + (i+0.5) * dx;
        double temp = temperature(u[i]);
        //std::cout<<"pressure is now "<<u[i][2]<<" at "<<i<<std::endl;
        out << x << " " << results[i][0] <<  " " << results[i][1] <<  " " << results[i][2] << " " << temp<<std::endl;

        // out<< r << results[i][0]<<results[i][1]<<results[i][2]<<"\n";
    }
    out.close();
}

void storeTimeData(const double pressure, int counter, std::array<double, 24> results){
    results[counter] = pressure ;
}

void PlotWithTime(std::array<double, 24> results){
    std::string filename = "SLIC.dat";
    std::ofstream out(filename);
    for(int i=0; i<24; i++){
        out<< i << " "<<results[i]<<"\n";
    }
}

//define a general flux function 

array flux_def(array x){
    array flux;
    array y = ConservativeToPrimative(x);
    double rho = x[0];
    double u = y[1];
    double p = y[2];

    flux[0] = x[1];               // mass: rho u
    flux[1] = rho * u * u + p;    // momentum: rho u² + p
    flux[2] = (x[2] + p) * u;     // energy: (E + p)u
    return flux;
}

//define a function to get the flux x is ui and y is ui+1 they are both arrays of length 3
//it will spit out an array of size 3 as well give each variable their flux

array getFlux(array  x , array y, double dx , double dt){
    //impliment general flux function 
    array  f_1 = flux_def(x); //f(ui)
    array  f_2 = flux_def(y); //f(i+1)
    array  uPlusHalf; //u_{i+1/2}

    for(int i=0; i<=2; ++i){
        uPlusHalf[i] = 0.5*(x[i] + y[i]) - 0.5*(dt/dx)*(f_2[i] - f_1[i]);
    }

    array RI_flux = flux_def(uPlusHalf); //richtmyer flux
    array LF_flux; //set up 3 length array for LF flux
    array FORCE_flux; //set up 3 length array for FORCE

    for (int i=0; i<=2; ++i) {
        LF_flux[i] = (0.5 * (dx/dt) * (x[i] - y[i])) + 0.5 * (f_1[i] + f_2[i]);
        FORCE_flux[i] = 0.5 * (LF_flux[i] + RI_flux[i]);
    }

    return FORCE_flux;
}

double computeTimeStep(const big_array& u , double C, double dx) {
    double maxSpeed = 0.0;

    for (const auto& state : u) {
        double rho = state[0];
        double mom = state[1];
        double E = state[2];

        array prim = ConservativeToPrimative(state);

        double u = prim[1];
        double pressure = prim[2];

        double sound_speed = SoundSpeed(state);
        double speed = std::abs(u) + sound_speed;

        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    double dt = C * dx / maxSpeed;
     return dt;
}

// Update source terms
big_array SourceTermUpdate(big_array u , double x0,double dx, double dt){
    big_array update(u.size());
    double alpha = 1.0;

    for(int i = 0; i < u.size(); i++) { 
    
        double r = x0 + (i+0.5) * dx;
        double rho = u[i][0];
        double mom = u[i][1];
        double E = u[i][2];

        array prim = ConservativeToPrimative(u[i]);
        double v = prim[1];
        double p = prim[2];

        // Source terms
        double S_rho = -alpha * rho * v / r;
        double S_mom = -alpha * rho * v * v / r;
        double S_E = -alpha * (E + p) * v / r;

        // Update using RK2 for source term only:

        // Stage 1
        array u_stage1;
        u_stage1[0] = rho + dt * S_rho;
        u_stage1[1] = mom + dt * S_mom;
        u_stage1[2] = E + dt * S_E;

        // Recompute primitives for stage 2
        array prim_stage1 = ConservativeToPrimative(u_stage1);
        double v1 = prim_stage1[1];
        double p1 = prim_stage1[2];

        // Stage 2 source terms
        double S_rho_2 = -alpha * u_stage1[0] * v1 / r;
        double S_mom_2 = -alpha * u_stage1[0] * v1 * v1 / r;
        double S_E_2 = -alpha * (u_stage1[2] + p1) * v1 / r;

        // Final update
        update[i][0] = rho + 0.5 * dt * (S_rho + S_rho_2);
        update[i][1] = mom + 0.5 * dt * (S_mom + S_mom_2);
        update[i][2] = E + 0.5 * dt * (S_E + S_E_2);
        
    }
    
    return update;
}

void applyBoundaryConditions(big_array& u){
    int n = u.size();
    // ----- REFLECTIVE ------
    u[0] = u[3];
    u[1] = u[2];
    // u[n - 1] = u[n - 4];
    // u[n - 2] = u[n - 3];
    u[0][1] = -u[3][1];
    u[1][1] = -u[2][1];
    // u[n - 1][1] = -u[n - 4][1];
    // u[n - 2][1] = -u[n - 3][1];

    // ----- TRANSMISSIVE -----
    // u[0] = u[3];
    // u[1] = u[2];
    u[n - 1] = u[n - 4];
    u[n - 2] = u[n - 3];

   
}

// a sub-diag, b main-diag, c sup-diag, d rhs, 
void ThomasAlorithm(const std::vector<double>& a, std::vector<double>& b,const std::vector<double>& c,  std::vector<double>& d) {
    int n = b.size();
    // Forward sweep
    for (int i = 1; i < n; ++i) {
        double m = a[i-1] / b[i-1];
        b[i] = b[i] - m * c[i-1];
        d[i] = d[i] - m * d[i-1];
    }
    // Back substitution
    d[n-1] = d[n-1] / b[n-1];
    for (int i = n - 2; i >= 0; --i) {
        d[i] = (d[i] - c[i] * d[i+1]) / b[i];
    }
}

//solve the poisson equation for magnetic potential
std::tuple<std::vector<double>, std::vector<double>> SolvePotential(double t,  double x0, double dx, double nCells, double alpha = 4137.95, double beta = 114866, double gamma = 10847100, double I0 = 218000){
    //set up initial data for current density with current I(t)
    double I = I0*(std::exp(-alpha * t) - std::exp(-beta * t))*(1 - std::exp(-gamma * t))*(1 - std::exp(-gamma * t));
    I = I0*std::exp(-alpha*t)*std::sin(beta*t);
    std::vector<double> J(nCells);
    for (int i =0; i<nCells; i++){
        double r = x0 + (i+0.5)*dx;
        J[i] = -I / (PI*r0*r0) * std::exp(-(r / r0)*(r / r0));
    }

    //coefficients of tridiagonal matrix
    std::vector<double> a(nCells - 1, -1.0);
    std::vector<double> b(nCells, 2.0);
    std::vector<double> c(nCells - 1, -1.0);
    std::vector<double> d(nCells);

    //set up RHS of the equation, d = -mu_0J*dx^2
    for (int i =0; i<nCells; i++){
        d[i] = -mu_0 * J[i] * dx * dx;
    }

    //add dirichlet boundary conditions (identity matrix at boundaries)
    b[0] = 1.0;
    c[0] = 0.0;
    d[0] = 0.0;
    b[nCells - 1] = 1.0;
    a[nCells - 2] = 0.0;
    d[nCells - 1] = 0.0;

    //solve the tridiagonal system using the thomas algorithm as in the paper
    ThomasAlorithm(a,b,c,d);

    std::vector<double> A(nCells);
    A = d;

    return {A,J};
}

//use finite difference to retrive the magnetic feild from it's potential using central differences
std::vector<double> MagneticFeild(std::vector<double> A, double dx){
    std::vector<double> B(A.size());

    //central differences
    for (int i=1; i<A.size()-1; i++){
        B[i] = - (A[i+1] - A[i-1]) / (2.0 * dx);
    }

    //endpoints
    B[0] = -(A[1] - A[0]) / dx;
    B[B.size() - 1] = -(A[A.size() - 1] - A[A.size() - 2]) / dx;

    return B;
}



//use euler method to update momentum
big_array momentumUpdate1(big_array u, double x0, double dx, double t, double dt){
    // std::vector<double> A(u.size());
    // std::vector<double> J(u.size());
    auto [A,J] = SolvePotential(t, x0, dx, u.size());
    std::vector<double> B(u.size());
    B = MagneticFeild(A,dx);

    big_array update(u.size());
    update = u;
    std::vector<double> cross(u.size());
    for (int i =0; i<cross.size(); i++){
        // J_z * B_{theta} (in -r direction)
        cross[i] = -J[i] * B[i]; 
    }

    int counter = 1000000;
    double dt_new = dt/counter;
    
    for(int i = 0; i < u.size(); i++) { 
    
        //update[i][1] = update[i][1] + dt*cross[i];
        
    
        // while(std::abs(1+cross[i]*dt_new) >= 1){ counter +=1; dt_new = dt/counter;std::cout<<counter<<std::endl;}
        int j = 0;
        while(j<=counter){
            update[i][1] = update[i][1] + dt_new*cross[i];
            j+=1;
        }
    

    }
    
    return update;
}
//use RK3 method to update momentum
big_array momentumUpdate(big_array u, double x0, double dx, double t, double dt) {
    auto [A, J] = SolvePotential(t, x0, dx, u.size());
    std::vector<double> B = MagneticFeild(A, dx);

    // Compute Lorentz force: cross[i] = -J[i] * B[i]
    std::vector<double> cross(u.size());
    for (int i = 0; i < u.size(); ++i) {
        cross[i] = -J[i] * B[i];
    }

    // Stage 1
    big_array u1 = u;
    for (int i = 0; i < u.size(); ++i) {
        u1[i][1] = u[i][1] + dt * cross[i];  // first Euler step
    }

    // Stage 2: Recompute force from u1 (assumes cross stays constant here, otherwise recompute)
    big_array u2 = u;
    for (int i = 0; i < u.size(); ++i) {
        u2[i][1] = 0.75 * u[i][1] + 0.25 * (u1[i][1] + dt * cross[i]);  // second stage
    }

    // Stage 3: Final RK3 step
    big_array update = u;
    for (int i = 0; i < u.size(); ++i) {
        update[i][1] = (1.0/3.0) * u[i][1] + (2.0/3.0) * (u2[i][1] + dt * cross[i]);  // final update
    }

    return update;
}
//use implicit method to update momentum
big_array momentumUpdate_Implicit(big_array u, double x0, double dx, double t, double dt, int max_iter = 10, double tol = 1e-6) {
    const int nCells = u.size();
    big_array update = u;

    // Initial guess: u^{n+1} ≈ u^n
    big_array guess = u;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Step 1: Compute J and B based on current guess
        auto [A, J] = SolvePotential(t + dt, x0, dx, nCells);  // using t + dt for implicitness
        std::vector<double> B = MagneticFeild(A, dx);

        // Step 2: Compute cross product force
        std::vector<double> cross(nCells);
        for (int i = 0; i < nCells; ++i) {
            cross[i] = -J[i] * B[i];
        }

        // Step 3: Update momentum using implicit formula
        double max_diff = 0.0;
        for (int i = 0; i < nCells; ++i) {
            double old_momentum = guess[i][1];
            double new_momentum = u[i][1] + dt * cross[i];

            guess[i][1] = new_momentum;

            double diff = std::abs(new_momentum - old_momentum);
            max_diff = std::max(max_diff, diff);
        }

        // Check for convergence
        // if (max_diff < tol) {
        //     std::cout << "[Implicit] Converged in " << iter + 1 << " iterations.\n";
        //     break;
        // }

        // if (iter == max_iter - 1) {
        //     std::cerr << "[Implicit] WARNING: Did not converge after " << max_iter << " iterations.\n";
        // }
    }

    update = guess;
    return update;
}



//use euler method to update energy
big_array energyUpdate1(big_array u, double x0, double dx, double t, double dt){
    auto [A,J] = SolvePotential(t, x0, dx, u.size());
    std::vector<double> B(u.size());
    B = MagneticFeild(A,dx);

    big_array update(u.size());
    update = u;
    std::vector<double> cross(u.size());
    for (int i =0; i<cross.size(); i++){
        // J_z * B_{theta} (in -r direction)
        cross[i] = -J[i] * B[i]; 
    }

    int counter = 1000000;
    double dt_new = dt/counter;

    for(int i = 0; i < u.size(); i++) { 
        double v = u[i][1] / u[i][0];

        
        // update[i][2] = update[i][2] + dt*v*cross[i];
        
    
        // while(std::abs(1+v*cross[i]*dt_new) >= 1){ counter +=1; dt_new = dt/counter;}
        int j = 0;
        while(j<=counter){
            update[i][2] = update[i][2] + dt_new*v*cross[i];
            j+=1;
        }
        
        
    }
    
    return update;
}
//use RK3 method to update energy
big_array energyUpdate(big_array u, double x0, double dx, double t, double dt) {
    auto [A, J] = SolvePotential(t, x0, dx, u.size());
    std::vector<double> B = MagneticFeild(A, dx);

    // Compute Lorentz force: cross[i] = -J[i] * B[i]
    std::vector<double> cross(u.size());
    for (int i = 0; i < u.size(); ++i) {
        cross[i] = -J[i] * B[i];
    }

    // Stage 1
    big_array u1 = u;
    for (int i = 0; i < u.size(); ++i) {
        u1[i][2] = u[i][2] + dt * cross[i]*u[i][1]/u[i][0];  // first Euler step
    }

    // Stage 2: Recompute force from u1 (assumes cross stays constant here, otherwise recompute)
    big_array u2 = u;
    for (int i = 0; i < u.size(); ++i) {
        u2[i][2] = 0.75 * u[i][2] + 0.25 * (u1[i][2] + dt * cross[i]*u1[i][1]/u1[i][0]);  // second stage
    }

    // Stage 3: Final RK3 step
    big_array update = u;
    for (int i = 0; i < u.size(); ++i) {
        update[i][2] = (1.0/3.0) * u[i][2] + (2.0/3.0) * (u2[i][2] + dt * cross[i]*u2[i][1]/u2[i][0]);  // final update
    }

    return update;
}
//use implicit metho to update energy
big_array energyUpdate_Implicit(big_array u, double x0, double dx, double t, double dt, int max_iter = 10, double tol = 1e-6) {
    const int nCells = u.size();
    big_array update = u;

    // Initial guess for energy: E^{n+1} ≈ E^n
    big_array guess = u;

    for (int iter = 0; iter < max_iter; ++iter) {
        // Step 1: Compute J and B from guess (at t + dt)
        auto [A, J] = SolvePotential(t + dt, x0, dx, nCells);
        std::vector<double> B = MagneticFeild(A, dx);

        // Step 2: Compute cross product force and energy source
        std::vector<double> energy_source(nCells);
        for (int i = 0; i < nCells; ++i) {
            double rho = guess[i][0];
            double momentum = guess[i][1];
            double velocity = momentum / rho;

            energy_source[i] = velocity * (-J[i] * B[i]);  // v · (J × B)
        }

        // Step 3: Update energy and check convergence
        double max_diff = 0.0;
        for (int i = 0; i < nCells; ++i) {
            double old_energy = guess[i][2];
            double new_energy = u[i][2] + dt * energy_source[i];

            guess[i][2] = new_energy;

            double diff = std::abs(new_energy - old_energy);
            max_diff = std::max(max_diff, diff);
        }

        // Check for convergence
        // if (max_diff < tol) {
        //     std::cout << "[Implicit Energy] Converged in " << iter + 1 << " iterations.\n";
        //     break;
        // }

        // if (iter == max_iter - 1) {
        //     std::cerr << "[Implicit Energy] WARNING: Did not converge after " << max_iter << " iterations.\n";
        // }
    }

    update = guess;
    return update;
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



//use euler to update energy with radiation
big_array thermalSourceTerm(big_array u,  double dt, double t, double dx, double x0, double nCells){
    //calculate the current density
    auto [A,J] = SolvePotential(t,x0,dx,nCells);

    //store the ellectrical conductivity at each time step
    std::vector<double> current_conductivity(u.size());

    //store the heat of formation at each space step
    std::vector<double> heat_of_formation(u.size());

    //store the temperature at each space step
    std::vector<double> current_temperature(u.size());

    //store the following velocity squared at each space step
    std::vector<double> current_velocity2(u.size());

    for(int i=0; i<u.size(); i++){
        //store the current electrical conductivity
        current_conductivity[i] = interpolate(u[i], electrical_conductivity);

        //array to store the current mass fractions at each space step
        std::array<double, 19> masses;
        for(int j =0; j<19; j++){
            masses[j]= interpolate(u[i],mass_fractions[j]);
            heat_of_formation[i] += masses[j] * heats_of_formation[j];
        }

        current_temperature[i] = temperature(u[i]);

        if(i != u.size() - 1){current_velocity2[i] = u[i+1][1]*u[i+1][1];}
        else{current_velocity2[i] = u[i][1] * u[i][1];}
    }

    big_array update(u.size());
    update = u;
    double kappa = 60;
    double stephen_botz = 5.67e-8;
    double radiation_term;

    for(int i=0; i<u.size(); i++){
        radiation_term = (stephen_botz*kappa*(std::pow(current_temperature[i], 4) - std::pow(T0,4)) - u[i][0]*(heat_of_formation[i] + 0.5*current_velocity2[i]));
        update[i][2] += dt * ((1.0/current_conductivity[i])*J[i]*J[i] - radiation_term);
    }
    
    return update;
    
}
//use implicit to update energy with radiation
big_array thermalSourceTerm_Implicit(big_array u, double dt, double t, double dx, double x0, double nCells, int max_iter = 10, double tol = 1e-6) {
    // Get current density
    auto [A, J] = SolvePotential(t + dt, x0, dx, nCells);  // use t + dt for implicit step

    big_array update = u;
    big_array guess = u;

    double kappa = 60;
    double stefan_boltzmann = 5.67e-8;

    for (int iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0;

        for (int i = 0; i < u.size(); ++i) {
            // --- Density and velocity ---
            double rho = u[i][0];
            double momentum = u[i][1];
            double velocity2 = momentum * momentum / (rho * rho);

            // --- Energy guess ---
            double E_guess = guess[i][2];
            double T_guess = temperature(guess[i]);

            // --- Interpolate conductivity ---
            double conductivity = interpolate(guess[i], electrical_conductivity);

            // --- Heat of formation ---
            double heat_of_formation = 0.0;
            for (int j = 0; j < 19; ++j) {
                double mass_frac = interpolate(guess[i], mass_fractions[j]);
                heat_of_formation += mass_frac * heats_of_formation[j];
            }

            // --- Radiation term ---
            double radiation = stefan_boltzmann * kappa * (std::pow(T_guess, 4) - std::pow(T0, 4));
            double loss_term = radiation - rho * (heat_of_formation + 0.5 * velocity2);

            // --- Joule heating ---
            double joule = (1.0 / conductivity) * J[i] * J[i];

            // --- Implicit update ---
            double E_new = u[i][2] + dt * (joule - loss_term);

            // --- Convergence check ---
            double diff = std::abs(E_new - E_guess);
            max_diff = std::max(max_diff, diff);

            guess[i][2] = E_new;  // update guess
        }

        if (max_diff < tol) {
            std::cout << "[Implicit Thermal] Converged in " << iter + 1 << " iterations.\n";
            break;
        }

        if (iter == max_iter - 1) {
            std::cerr << "[Implicit Thermal] WARNING: Did not converge after " << max_iter << " iterations.\n";
        }
    }

    update = guess;
    return update;
}
//use implicit to update energy with radiation
big_array thermalSourceTerm_Newton(big_array u, double dt, double t, double dx, double x0, double nCells, int max_iter = 20, double tol = 1e-6) {
    auto [A, J] = SolvePotential(t + dt, x0, dx, nCells);  // use t + dt for implicit
    big_array update = u;

    double kappa = 60;
    double stefan_boltzmann = 5.67e-8;
    double epsilon = 1e-5;  // for finite difference

    for (int i = 0; i < u.size(); ++i) {
        double E_old = u[i][2];  // E^n
        double E = E_old;        // initial guess E^0

        double rho = u[i][0];
        double momentum = u[i][1];
        double velocity2 = (momentum / rho) * (momentum / rho);

        // Precompute J²
        double J2 = J[i] * J[i];

        // Compute constant mass-fraction-based heat of formation
        double heat_of_formation = 0.0;
        for (int j = 0; j < 19; ++j) {
            double mass_frac = interpolate(u[i], mass_fractions[j]);
            heat_of_formation += mass_frac * heats_of_formation[j];
        }

        int iter = 0;
        double diff = 1e9;

        while (iter < max_iter && diff > tol) {
            // T from energy guess
            double T = temperature({rho, momentum, E});

            // σ from E guess
            double sigma = interpolate({rho, momentum, E}, electrical_conductivity);

            // Radiation loss term S_T
            double S_T = stefan_boltzmann * kappa * (std::pow(T, 4) - std::pow(T0, 4)) - rho * (heat_of_formation + 0.5 * velocity2);

            // f(E)
            double f = E - E_old - dt * ((1.0 / sigma) * J2 - S_T);

            // f'(E) via finite difference
            double E_eps = E + epsilon;

            // T_eps from E+ε
            double T_eps = temperature({rho, momentum, E_eps});
            double sigma_eps = interpolate({rho, momentum, E_eps}, electrical_conductivity);
            double S_T_eps = stefan_boltzmann * kappa * (std::pow(T_eps, 4) - std::pow(T0, 4)) - rho * (heat_of_formation + 0.5 * velocity2);
            double f_eps = E_eps - E_old - dt * ((1.0 / sigma_eps) * J2 - S_T_eps);

            double df_dE = (f_eps - f) / epsilon;

            // Newton step
            double E_new = E - f / df_dE;
            diff = std::abs(E_new - E);
            E = E_new;
            ++iter;
        }


        update[i][2] = E;  // set updated energy
    }

    return update;
}


int main() { 
    int nCells = 100; //the distance between points is 0.01
    double x0 = 0.0;
    double x1 = 0.2;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 1.5e-4;
    double C = 0.8;
    double omega = 0;

    // Allocate matrices with 2 extra points for transmissive BCs
    big_array  u(nCells+4);
    big_array  flux(u.size());
    big_array  uBarL(u.size());
    big_array  uBarR(u.size());
    big_array  uBarHalfL(u.size());
    big_array  uBarHalfR(u.size());
    big_array  uPlus1(u.size());
    double dx = (x1 - x0) / nCells; //the space steps 
    double time;

    // Initial conditions!

    for(int i = 0; i < u.size(); i++) {
        // x 0 is at point i=1/2
        double x = x0 + (i+0.5) * dx;
        std::array<double, 3> prim;
        // if(x <= 0.4) {
        //     prim[0] = 1; // Density
        //     prim[1] = 0*std::pow(10,2.5); // Velocity
        //     prim[2] = 1*std::pow(10,5); // Pressure
        //     } else {
        //     prim[0] = 0.125; // Density
        //     prim[1] = 0*std::pow(10,2.5); // Velocity
        //     prim[2] = 0.1*std::pow(10,5); // Pressure
        // }

        prim[0] = 1.225; // Density
        prim[1] = 0*std::pow(10,2.5); // Velocity
        prim[2] = 101325 + 2e6*std::exp(-(x/r0)*(x/r0)); // Pressure

        u[i] = PrimativeToConservative(prim);
        // array v = ConservativeToPrimative(u[i]);
    }

    double dt = computeTimeStep(u , C , dx); //the time steps
    double t = tStart;
    int counter =0;
    do {
        // Compute the stable time step for this iteration

        dt = computeTimeStep(u , C , dx); 
        t = t + dt;
        std::cout<<"t= "<<t<<" dt= "<<dt<<std::endl;

        // Update source terms
        u = SourceTermUpdate(u,x0,dx,dt);
        //update resistive source terms
        u = momentumUpdate_Implicit(u, x0, dx, t, 0.5*dt);
        u = energyUpdate_Implicit(u, x0, dx, t, 0.5*dt);

        // Apply boundary conditions
        applyBoundaryConditions(u);

        //find ubar


        for(int i=1; i<=nCells+2; ++i){
            for(int j=0; j<=2; ++j){
                double DeltaPlus = u[i+1][j] - u[i][j];
                double DeltaMinus = u[i][j] - u[i-1][j];
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

                uBarL[i][j] = u[i][j] - 0.5 * xi * Delta;
                uBarR[i][j] = u[i][j] + 0.5 * xi * Delta;
                
                
            }
        }

        for(int i=1; i<=nCells+2; ++i){
            for(int j=0; j<=2; ++j){
                uBarHalfL[i][j] = uBarL[i][j] - 0.5*(dt/dx)*(flux_def(uBarR[i] )[j]-flux_def(uBarL[i] )[j]);
                uBarHalfR[i][j] = uBarR[i][j] - 0.5*(dt/dx)*(flux_def(uBarR[i] )[j]-flux_def(uBarL[i] )[j]);
                
            }
        }

        
        // Apply boundary conditions
        applyBoundaryConditions(uBarHalfL);
        applyBoundaryConditions(uBarHalfR);


        for(int i = 0; i <= nCells+3; i++) { //Define the fluxes
            // flux[i] corresponds to cell i+1/2 
            flux[i] = getFlux( uBarHalfR[i], uBarHalfL[i+1] , dx , dt);
        }

        //the below has a i-1 flux which means we need to define a flux at 0 so make sure the above ^ starts at 0! this is because we have another edge with the number of cells (like the walls)

        for(int i = 1; i <= nCells+3; i++) { //Update the data
            for(int j=0; j<=2; ++j){
                uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
            }
        }
    
        // Now replace u with the updated data for the next time step
        u = uPlus1;
        

        // Output data at specific time steps
        while (t >= 1e-5 * counter) {
            big_array results(u.size());

            // for (int i = 0; i <= results.size() - 1; ++i) {
            //     results[i] = ConservativeToPrimative(u[i]);
            // }

            SaveUnwrappedData(u, x0, dx, nCells, counter);
            // storeTimeData(u[u.size()-1][2], counter, results);
            std::cout << "Saved frame: " << counter << std::endl;
            counter += 1;
        }

        //radiation source term
        u = thermalSourceTerm_Newton(u, dt,t,dx,x0,nCells);
        applyBoundaryConditions(u);

        //update resistive source terms
        u = momentumUpdate_Implicit(u, x0, dx, t, 0.5*dt);
        u = energyUpdate_Implicit(u, x0, dx, t, 0.5*dt);
        
                
               
    } while (t < tStop);

    applyBoundaryConditions(u);

    // PlotWithTime(results);
    //still need to convert it back to primitive

    //define final results

    big_array results(u.size());

    for(int i=0; i<= results.size() -1; ++i){
        results[i] = ConservativeToPrimative(u[i]);
    }


    //output
    std::string filename = "euler.dat";
    std::ofstream output(filename);
    for (int i = 2; i <= nCells+1; ++i) {
        double x = x0 + (i+0.5) * dx;
        double temp = temperature(u[i]);
        //std::cout<<"pressure is now "<<u[i][2]<<" at "<<i<<std::endl;
        output << x << " " << results[i][0] <<  " " << results[i][1] <<  " " << results[i][2] << " " << temp<<std::endl;
    }

    //wrap data around the r=0 axis
    SaveUnwrappedData(u, x0, dx, nCells, counter);
    
}