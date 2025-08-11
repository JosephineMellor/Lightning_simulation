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
typedef std::array<double,8> array;
typedef std::vector<std::array<double,8>> big_array;

const double PI = 3.141592653589793;

//function to read in the data from the tabulated equation of state
std::tuple<data_vec , data_vec , data_table , data_table , data_table , data_table , data_table> Plasma19(){
    
    //load the file
    std::ifstream file("Plasma19_EoS.txt"); 

    //set up the lines, starting row, and data arrays
    std::string line;
    int row = 0;
    data_vec densities;
    data_vec pressures;

    //we want to save sound speeds, energies, temperatures, electrical conductivity and thermal conductivity in tables
    data_table sound_speeds(500, std::vector<double>(500));
    data_table energies(500, std::vector<double>(500));
    data_table temperatures(500, std::vector<double>(500));
    data_table electrical_conductivity(500, std::vector<double>(500));
    data_table thermal_conductivity(500, std::vector<double>(500));
    
    //now loop over the data file row by row, we only want the fifth and sixth row for densities and pressures
    while(std::getline(file , line)){
        row++;

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

        //thermal conductivity
        if( row >= 2018 && row< 2518){
            std::istringstream iss(line); 

            big_vector row_data; 
            double val;
            while (iss >> val) {
                row_data.push_back(val);
            }

            thermal_conductivity[row-2018] = std::move(row_data);
        }
    }


    return{densities , pressures , sound_speeds , energies , temperatures , electrical_conductivity , thermal_conductivity};
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

//update source terms
// std::array<double, 3> MomentumUpdate(std::array<double , 8>u){
//     double dt = 0.001;
//     std::array<double , 8>uBar;
//     int count;
//     while(count < 100){
//         for(int i=0; i<u.size() ; i++){
//             uBar[i] = u[i] + dt ;
//         }
//     }
// }
std::array<double , 8> PrimitiveToConservative(const std::array<double , 8>& u ){
    //start with variables that dont require equation of state
    std::array<double ,8> v;
    v[0] = u[0]; //density
    v[1] = u[0] * u[1]; //vx
    v[2] = u[0] * u[2]; //vy
    v[3] = u[0] * u[3]; //vz
    v[5] = u[5]; //Bx
    v[6] = u[6]; //By
    v[7] = u[7]; //Bz

    double B2 = (u[5]*u[5] + u[6]*u[6] + u[7]*u[7]);

    //find index for density with binary search
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
    if(density_bottom == density_top){density_ratio =1;}
    

    //find index for pressure also with binary search
    point1 =0;
    point2 =499;
    half;
    int pressure_index;
    double p = u[4];
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

    v[4] = energy*rho + 0.5*rho*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]) + 0.5 * B2;
    v[4] = u[4] / (2.0-1) + 0.5*u[0]*(u[1]*u[1] + u[2]*u[2] + u[3]*u[3]) + 0.5*(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]); //energy
    return v;
}
std::array<double, 8> ConservativeToPrimitive(const std::array<double , 8>& u){
    //start with variables that dont require equation of state
    std::array<double , 8> v;
    v[0] = u[0];
    v[1] = u[1] / u[0];
    v[2] = u[2] / u[0];
    v[3] = u[3] / u[0];
    v[5] = u[5];
    v[6] = u[6];
    v[7] = u[7];

    double B2 = (u[5]*u[5] + u[6]*u[6] + u[7]*u[7]);

    //find our index for density with binary search
    int point1 =0;
    int point2 =499;
    int half;
    int density_index;
    double rho = v[0];
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
    double e = (u[3] - 0.5*rho*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]) - 0.5*B2) / rho; //internal energy

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

    v[4] = p;
    v[4] = (u[4] - 0.5*v[0]*(v[1]*v[1] + v[2]*v[2] + v[3]*v[3]) - 0.5*(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]))*(2.0 -1);//pressure

    return v;
}

double SoundSpeed(const std::array<double,8> u){
    std::array<double , 8> prim = ConservativeToPrimitive(u);

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
    sound_speed = std::sqrt(2.0 * prim[4] / rho);
    
    return sound_speed;
}

std::array<double , 8> FluxDef(const std::array<double , 8>& u ){//conservative variables go into this function
    std::array<double , 8> v = ConservativeToPrimitive(u);
    std::array<double , 8> flux;
    flux[0] = u[1]; //rho v_x
    flux[1] = u[0]*v[1]*v[1] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ) - u[5]*u[5]; 
    flux[2] = u[0]*v[1]*v[2] - u[5]*u[6];
    flux[3] = u[0]*v[1]*v[3] - u[5]*u[7];
    flux[4] = (u[4] + v[4] + 0.5*( u[5]*u[5] + u[6]*u[6] + u[7]*u[7] ))*v[1] - (v[1]*v[5]+v[2]*v[6]+v[3]*v[7])*v[5];
    flux[5] = 0;
    flux[6] = v[6]*v[1]-v[5]*v[2];
    flux[7] = v[7]*v[1]-v[5]*v[3];
    return flux;
}
double ComputeTimeStep(const std::vector<std::array<double,8>>& u , double C, double dx) {
    double maxSpeed = 0.0;

    for (const auto& state : u) {
        double rho = state[0];
        double mom_x = state[1];
        double mom_y = state[2];
        double mom_z = state[3];
        double E = state[4];
        double Bx = state[5];
        double By = state[6];
        double Bz = state[7];

        double u_x = mom_x / rho;
        double u_y = mom_y / rho;
        double u_z = mom_z / rho;
        double intermediate = 0.5 * rho *( u_x * u_x + u_y*u_y + u_z*u_z) + 0.5*(Bx*Bx + By*By + Bz*Bz);
        double BmagSquared = Bx*Bx + By*By + Bz*Bz;
        double pressure = ConservativeToPrimitive(state)[4];

        double sound_speed = SoundSpeed(state);
        double alfven_speed = std::abs(Bx) / std::sqrt(rho);
        double slow_ma_speed = std::sqrt( 0.5*(sound_speed*sound_speed + (BmagSquared/rho) - std::sqrt((sound_speed*sound_speed + BmagSquared/rho)*(sound_speed*sound_speed + BmagSquared/rho) - 4.0*sound_speed*sound_speed*Bx*Bx / rho)));
        
        double fast_ma_speed = std::sqrt( 0.5*(sound_speed*sound_speed + (BmagSquared/rho) + std::sqrt((sound_speed*sound_speed + BmagSquared/rho)*(sound_speed*sound_speed + BmagSquared/rho) - 4.0*sound_speed*sound_speed*Bx*Bx / rho)));
        
        double speed = std::abs(u_x) + fast_ma_speed;

        if (speed > maxSpeed) {
            maxSpeed = speed;
        }
    }

    double dt = C * dx / maxSpeed;
     return dt;
}
//define a function to get the flux x is ui and y is ui+1 they are both arrays of length 3
//it will spit out an array of size 3 as well give each variable their flux

std::array<double , 8> uHLL ( std::array<double , 8>x , std::array<double , 8> y){
        //define all the variables
    double rhoL = x[0];
    double mom_xL = x[1];
    double mom_yL = x[2];
    double mom_zL = x[3];
    double EL = x[4];
    double BxL = x[5];
    double ByL = x[6];
    double BzL = x[7];

    double u_xL = mom_xL / rhoL;
    double u_yL = mom_yL / rhoL;
    double u_zL = mom_zL / rhoL;
    double intermediateL = 0.5 * rhoL *( u_xL * u_xL + u_yL*u_yL + u_zL*u_zL) + 0.5*(BxL*BxL + ByL*ByL + BzL*BzL);
    double BmagSquaredL = BxL*BxL + ByL*ByL + BzL*BzL;
    double pressureL = ConservativeToPrimitive(x)[4];
    double pTL = pressureL + 0.5*BmagSquaredL;

    double sound_speedL = SoundSpeed(x);
    double alfven_speedL = std::abs(BxL) / std::sqrt(rhoL);
    double slow_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) - std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    
    double rhoR = y[0];
    double mom_xR = y[1];
    double mom_yR = y[2];
    double mom_zR = y[3];
    double ER = y[4];
    double BxR = y[5];
    double ByR = y[6];
    double BzR = y[7];

    double u_xR = mom_xR / rhoR;
    double u_yR = mom_yR / rhoR;
    double u_zR = mom_zR / rhoR;
    double intermediateR = 0.5 * rhoR * (u_xR * u_xR + u_yR * u_yR + u_zR * u_zR) + 0.5 * (BxR * BxR + ByR * ByR + BzR * BzR);
    double BmagSquaredR = BxR * BxR + ByR * ByR + BzR * BzR;
    double pressureR = ConservativeToPrimitive(y)[4];
    double pTR = pressureR + 0.5*BmagSquaredR;

    double sound_speedR = SoundSpeed(y);
    double alfven_speedR = std::abs(BxR) / std::sqrt(rhoR);
    double slow_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) - std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));

    double sL = std::min(u_xL , u_xR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_xL , u_xR) + std::max(fast_ma_speedL , fast_ma_speedR);

    std::array<double , 8> uStar;
    auto FluxX = FluxDef(x);
    auto FluxY = FluxDef(y);
    for (int i = 0; i < 8; ++i) {
        uStar[i] = (sR * y[i] - sL * x[i] - FluxY[i] + FluxX[i]) / (sR - sL);
    }

    std::array<double , 8> uHLL;
    if (sL >= 0) {
        uHLL = x;
    } else if (sL < 0 && sR > 0) {
        uHLL = uStar;
    } else { // sR <= 0
        uHLL = y;
    }

    return uHLL;
}

//using the minbee limiter


std::array<double , 8> getFlux ( std::array<double , 8>x , std::array<double , 8> y ){
    //define all the variables
    double rhoL = x[0];
    double mom_xL = x[1];
    double mom_yL = x[2];
    double mom_zL = x[3];
    double EL = x[4];
    double BxL = x[5];
    double ByL = x[6];
    double BzL = x[7];

    double u_xL = mom_xL / rhoL;
    double u_yL = mom_yL / rhoL;
    double u_zL = mom_zL / rhoL;
    double intermediateL = 0.5 * rhoL *( u_xL * u_xL + u_yL*u_yL + u_zL*u_zL) + 0.5*(BxL*BxL + ByL*ByL + BzL*BzL);
    double BmagSquaredL = BxL*BxL + ByL*ByL + BzL*BzL;
    double pressureL = ConservativeToPrimitive(x)[4];
    double pTL = pressureL + 0.5*BmagSquaredL;

    double sound_speedL = SoundSpeed(x);
    double alfven_speedL = std::abs(BxL) / std::sqrt(rhoL);
    double slow_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) - std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    double fast_ma_speedL = std::sqrt( 0.5*(sound_speedL*sound_speedL + (BmagSquaredL/rhoL) + std::sqrt((sound_speedL*sound_speedL + BmagSquaredL/rhoL)*(sound_speedL*sound_speedL + BmagSquaredL/rhoL) - 4.0*sound_speedL*sound_speedL*BxL*BxL / rhoL)));
    
    double rhoR = y[0];
    double mom_xR = y[1];
    double mom_yR = y[2];
    double mom_zR = y[3];
    double ER = y[4];
    double BxR = y[5];
    double ByR = y[6];
    double BzR = y[7];

    double u_xR = mom_xR / rhoR;
    double u_yR = mom_yR / rhoR;
    double u_zR = mom_zR / rhoR;
    double intermediateR = 0.5 * rhoR * (u_xR * u_xR + u_yR * u_yR + u_zR * u_zR) + 0.5 * (BxR * BxR + ByR * ByR + BzR * BzR);
    double BmagSquaredR = BxR * BxR + ByR * ByR + BzR * BzR;
    double pressureR = ConservativeToPrimitive(y)[4];
    double pTR = pressureR + 0.5*BmagSquaredR;

    double sound_speedR = SoundSpeed(y);
    double alfven_speedR = std::abs(BxR) / std::sqrt(rhoR);
    double slow_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) - std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));
    double fast_ma_speedR = std::sqrt(0.5 * (sound_speedR * sound_speedR + (BmagSquaredR / rhoR) + std::sqrt((sound_speedR * sound_speedR + BmagSquaredR / rhoR) * (sound_speedR * sound_speedR + BmagSquaredR / rhoR) - 4.0 * sound_speedR * sound_speedR * BxR * BxR / rhoR)));

    double sL = std::min(u_xL , u_xR) - std::max(fast_ma_speedL , fast_ma_speedR);
    double sR = std::max(u_xL , u_xR) + std::max(fast_ma_speedL , fast_ma_speedR);

    //double qStar = (rhoR*u_xR*(sR - u_xR) - rhoL*u_xL*(sL - u_xL) + pTL - pTR - BxL*BxL + BxR*BxR) / ( rhoR*(sR - u_xR) - rhoL*(sL - u_xL));
    double qStar = (pTR - pTL + rhoL * u_xL * (sL - u_xL) - rhoR * u_xR * (sR - u_xR)) / (rhoL * (sL - u_xL) - rhoR * (sR - u_xR));

    std::array<double , 8> u_HLL = uHLL(x,y);
    double BxHLL = u_HLL[5];
    double ByHLL = u_HLL[6];
    double BzHLL = u_HLL[7];

    double rhoHLL = u_HLL[0];
    double u_xHLL = u_HLL[1] / rhoHLL;  
    double u_yHLL = u_HLL[2] / rhoHLL;  
    double u_zHLL = u_HLL[3] / rhoHLL;  

    double BxStarL = BxL;
    double ByStarL = ByHLL;
    double BzStarL = BzHLL;

    double pStarL = rhoL*(sL - u_xL)*(qStar - u_xL) + pTL - BxL*BxL + BxStarL*BxStarL;
    //double pStarL = rhoL*(sL - u_xL)*(qStar - u_xL) + pressureL - 0.5*BxL*BxL + 0.5*BxStarL*BxStarL;

    double rhoStarL = rhoL*(sL - u_xL)/(sL - qStar);
    double mom_xStarL = rhoStarL * qStar;
    double mom_yStarL = mom_yL*(sL-u_xL)/(sL - qStar) - (BxStarL*ByStarL - BxL*ByL)/(sL - qStar);
    double mom_zStarL = mom_zL*(sL-u_xL)/(sL - qStar) - (BxStarL*BzStarL - BxL*BzL)/(sL - qStar);
    double u_xStarL = mom_xStarL / rhoStarL;
    double u_yStarL = mom_yStarL / rhoStarL;
    double u_zStarL = mom_zStarL / rhoStarL;

    double vDotB_StarL = BxStarL * u_xStarL + ByStarL * u_yStarL + BzStarL * u_zStarL;
    double vDotB_L = BxL * u_xL + ByL * u_yL + BzL * u_zL;
    double energyStarL = (EL * (sL - u_xL) - pTL * u_xL + pStarL * qStar + BxL * (vDotB_L - vDotB_StarL)) / (sL - qStar);


    std::array<double ,8> uStarL;
    uStarL[0] = rhoStarL;
    uStarL[1] = mom_xStarL;
    uStarL[2] = mom_yStarL;
    uStarL[3] = mom_zStarL;
    uStarL[4] = energyStarL;
    uStarL[5] = BxStarL;
    uStarL[6] = ByStarL;
    uStarL[7] = BzStarL;


    double BxStarR = BxR;
    double ByStarR = ByHLL;
    double BzStarR = BzHLL;

    //double pStarR = rhoR*(sR - u_xR)*(qStar - u_xR) + pressureR - BxR*BxR + BxStarR*BxStarR;
    double pStarR = rhoR*(sR - u_xR)*(qStar - u_xR) + pTR - 0.5*BxR*BxR + 0.5*BxStarR*BxStarR;

    double rhoStarR = rhoR*(sR - u_xR)/(sR - qStar);
    double mom_xStarR = rhoStarR * qStar;
    double mom_yStarR = mom_yR*(sR-u_xR)/(sR - qStar) - (BxStarR*ByStarR - BxR*ByR)/(sR - qStar);
    double mom_zStarR = mom_zR*(sR-u_xR)/(sR - qStar) - (BxStarR*BzStarR - BxR*BzR)/(sR - qStar);
    double u_xStarR = mom_xStarR / rhoStarR;
    double u_yStarR = mom_yStarR / rhoStarR;
    double u_zStarR = mom_zStarR / rhoStarR;

    double vDotB_StarR = BxStarR * u_xStarR + ByStarR * u_yStarR + BzStarR * u_zStarR;
    double vDotB_R = BxR * u_xR + ByR * u_yR + BzR * u_zR;
    double energyStarR = (ER * (sR - u_xR) - pTR * u_xR + pStarR * qStar + BxR * (vDotB_R - vDotB_StarR)) / (sR - qStar);


    std::array<double ,8> uStarR;
    uStarR[0] = rhoStarR;
    uStarR[1] = mom_xStarR;
    uStarR[2] = mom_yStarR;
    uStarR[3] = mom_zStarR;
    uStarR[4] = energyStarR;
    uStarR[5] = BxStarR;
    uStarR[6] = ByStarR;
    uStarR[7] = BzStarR;

    std::array<double , 8> FluxL = FluxDef(x);
    std::array<double , 8> FluxR = FluxDef(y);
    std::array<double , 8> FluxStarL = FluxDef(uStarL);
    std::array<double , 8> FluxStarR = FluxDef(uStarR);
    std::array<double , 8> FluxHLLC;

    if (sL >= 0) {
        FluxHLLC = FluxL;
    } else if (sL <= 0 && qStar >= 0) {
        for (int i = 0; i < 8; ++i)
            FluxHLLC[i] = FluxL[i] + sL * (uStarL[i] - x[i]);
    } else if (qStar <= 0 && sR >= 0) {
        for (int i = 0; i < 8; ++i)
            FluxHLLC[i] = FluxR[i] + sR * (uStarR[i] - y[i]);
    } else {
        FluxHLLC = FluxR;
    }
    return FluxHLLC;
}


void applyBoundaryConditions(big_array& u){
    int n = u.size();
    // ----- REFLECTIVE ------
    // u[0] = u[3];
    // u[1] = u[2];
    // u[n - 1] = u[n - 4];
    // u[n - 2] = u[n - 3];
    // u[0][1] = -u[3][1];
    // u[1][1] = -u[2][1];
    // u[n - 1][1] = -u[n - 4][1];
    // u[n - 2][1] = -u[n - 3][1];

    // ----- TRANSMISSIVE -----
    u[0] = u[3];
    u[1] = u[2];
    u[n - 1] = u[n - 4];
    u[n - 2] = u[n - 3];

}



int main() { 
    int nCells = 800; //the distance between points is 0.01
    double x0 = 0.0;
    double x1 = 800;
    double tStart = 0.0; //set the start and finish time steps the same
    double tStop = 80;
    double C = 1.0;
    double omega = 0;

    // Allocate matrices with 2 extra points for transmissive BCs
    std::vector<std::array<double,8>> u(nCells+4);
    std::vector<std::array<double,8>> uPlus1(u.size());
    std::vector<std::array<double,8>> flux(u.size());
    std::vector<std::array<double,8>> uBarL(u.size());
    std::vector<std::array<double,8>> uBarR(u.size());
    std::vector<std::array<double,8>> uBarHalfL(u.size());
    std::vector<std::array<double,8>> uBarHalfR(u.size());
    double dx = (x1 - x0) / nCells; //the space steps 

    // Initial conditions!

    for(int i = 0; i < u.size(); i++) {
        // x 0 is at point i=1/2
        double x = x0 + (i-0.5) * dx;
        std::array<double, 8> prim;
        if(x <= 400) {
            prim[0] = 1; // Density
            prim[1] = 0; // Velocity
            prim[2] = 0;
            prim[3] = 0;
            prim[4] = 1; // pressure
            prim[5] = 0.75; // magnetic field
            prim[6] = 1;
            prim[7] = 0; 
            } else {
            prim[0] = 0.125; // Density
            prim[1] = 0; // Velocity
            prim[2] = 0;
            prim[3] = 0;
            prim[4] = 0.1; // pressure
            prim[5] = 0.75; // magnetic field
            prim[6] = -1;
            prim[7] = 0; 
        }

        u[i] = PrimitiveToConservative(prim);
    }

    

    double dt = ComputeTimeStep(u , C , dx); //the time steps

    double t = tStart;
    do {

        dt = ComputeTimeStep(u , C , dx); 
        t = t + dt;
        std::cout<<"t= "<< t<< " dt= "<< dt<< std::endl;

        applyBoundaryConditions(u);

        //find ubar


        for(int i=1; i<=u.size() - 2; ++i){
            for(int j=0; j<8; ++j){
                double DeltaPlus = u[i+1][j] - u[i][j];
                double DeltaMinus = u[i][j] - u[i-1][j];
                // std::cout<< DeltaMinus << " " << DeltaPlus << std::endl;
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

        for(int i=1; i<=u.size() - 2; ++i){
            for(int j=0; j<8; ++j){
                uBarHalfL[i][j] = uBarL[i][j] - 0.5*(dt/dx)*(FluxDef(uBarR[i])[j]-FluxDef(uBarL[i])[j]);
                uBarHalfR[i][j] = uBarR[i][j] - 0.5*(dt/dx)*(FluxDef(uBarR[i])[j]-FluxDef(uBarL[i])[j]);
                
            }
        }
        
        applyBoundaryConditions(uBarHalfL);
        applyBoundaryConditions(uBarHalfR);

        for(int i = 0; i < u.size() - 1; i++) { //Define the fluxes
            // flux[i] corresponds to cell i+1/2 
            flux[i] = getFlux( uBarHalfR[i], uBarHalfL[i+1]);
        }

        //the below has a i-1 flux which means we need to define a flux at 0 so make sure the above ^ starts at 0! this is because we have another edge with the number of cells (like the walls)

        for(int i = 1; i <= u.size() - 1; i++) { //Update the data
            for(int j=0; j<8; ++j){
                uPlus1[i][j] = u[i][j] - (dt/dx) * (flux[i][j] - flux[i-1][j]);
            }
        }
    
        // Now replace u with the updated data for the next time step

        
        u = uPlus1;
    } while (t < tStop);

    //still need to convert it back to primitive
    //define final results

    std::vector<std::array<double,8>> results(u.size());

    for(int i=0; i<= results.size() -1; ++i){
        results[i] = ConservativeToPrimitive(u[i]);
    }


    // Output the results
    std::ofstream output("MHD.dat");
    for (int i = 0; i <= u.size() - 1; ++i) {
        double x = x0 + (i - 1) * dx;
        output << x << " " << results[i][0] <<  " " << results[i][1] <<  " " << results[i][2] << " " << results[i][3] <<" " << results[i][4] <<" " << results[i][5] <<" " << results[i][6] <<" " << results[i][7] << std::endl;
        // std::cout << x << " " << u[i][0] <<  " " << u[i][1] <<  " " << u[i][2] << std::endl;
    }
}