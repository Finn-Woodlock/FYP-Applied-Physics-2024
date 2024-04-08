#include <Eigen>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>


namespace Tools
{

    typedef Eigen::Matrix<long double,Eigen::Dynamic,Eigen::Dynamic> DataMat;
    typedef Eigen::Matrix<std::string,Eigen::Dynamic,Eigen::Dynamic> DataMatString;

    typedef Eigen::Matrix<long double,1,Eigen::Dynamic> DataRowVector;
    typedef Eigen::Matrix<std::string,1,Eigen::Dynamic> DataRowVectorString;

        //----------------------------------------------------------------------------------------------------------

    void WriteMatrix(std::string dir,std::string file, std::string FirstRow = "", DataMat data = DataMat(0,0), DataMatString strdata = DataMatString(0,0) )
    {
        
        std::ofstream outdata(dir+file);

        if ( outdata.is_open() ) 
        { 
            outdata << FirstRow << std::endl;

            int i;
            int j;
            for (i=0; i<data.cols(); ++i)
            {
                for (j=0; j<data.rows(); ++j)
                {
                    if (j > 0) {outdata << " , ";}
                    outdata << data(j,i);
                }
                
                for (j=data.rows(); j<data.rows()+strdata.rows(); ++j)
                {
                    if (j > 0) {outdata << " , ";}
                    outdata << strdata(j-data.rows(),i);
                }
                outdata << std::endl;
            }
            outdata.close();
        }
        else
        {
            std::cout << "Error: file could not be opened" << std::endl;
            std::cerr << "Error: file could not be opened" << std::endl;
            exit(1);
        }

        return;
    }


        //----------------------------------------------------------------------------------------------------------

    std::vector<double> linspace(double start_in, double end_in, int num_in)
    {

    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1) 
        {
        linspaced.push_back(start);
        return linspaced;
        }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
        {
        linspaced.push_back(start + delta * i);
        }
    linspaced.push_back(end); 
    return linspaced;
    }


        //----------------------------------------------------------------------------------------------------------

    Eigen::Matrix2Xd PointArray2D(int x,int y){
        Eigen::Matrix2Xd p(2, x * y);
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)
            {
                p(0, j * x + i) = i;
                p(1, j * x + i) = j;
            }
        }
        return p;
    }

        //----------------------------------------------------------------------------------------------------------

    Eigen::Matrix3Xd PointArray3D(int x, int y, int z){
        Eigen::Matrix3Xd p(3, x * y * z);
        for (int i = 0; i < x; i++)
        {
            for (int j = 0; j < y; j++)
            {
                for (int k = 0; k < z; k++)
                {
                    p(0, k * x * y + j * x + k) = i;
                    p(1, k * x * y + j * x + k) = j;
                    p(2, k * x * y + j * x + k) = k;
                }
            }
        }
        return p;
    }

        //----------------------------------------------------------------------------------------------------------

    Eigen::Matrix2Xd PointArray2DRand(int x, double min, double max){
        Eigen::Matrix2Xd p(2, x);
        for (int i = 0; i < x; i++)
        {
            Eigen::VectorXd m = Eigen::VectorXd::Random(2,1); 
            m = (m + Eigen::VectorXd::Constant(2,1.))*(max-min)/2.;
            m = (m + Eigen::VectorXd::Constant(2,min));
            p.col(i) = m;
        }
        return p;
    }

        //----------------------------------------------------------------------------------------------------------

    Eigen::Matrix3Xd PointArray3DRand(int x, double min, double max){
        Eigen::Matrix3Xd p(3, x);
        for (int i = 0; i < x; i++)
        {
            Eigen::VectorXd m = Eigen::VectorXd::Random(3,1);
            m = (m + Eigen::VectorXd::Constant(3,1.))*(max-min)/2.;
            m = (m + Eigen::VectorXd::Constant(3,min));
            p.col(i) = m;
        }
        return p;
    }
   
    class BravaisLattice
    {
    public:
        BravaisLattice()
        {
            return;
        }

        static Eigen::Matrix2d tP2D()
        {
            Eigen::Matrix2d out;
            out << 1,0,
                   0,1;
            return out;
        }

        static Eigen::Matrix2d tC2D()
        {
            Eigen::Matrix2d out;
            out << 1,0.5,
                   0,0.5;
            return out;
        }

        static Eigen::Matrix2d hP2D()
        {
            Eigen::Matrix2d out;
            out << 1,std::cos(M_PI/3),
                   0,std::sin(M_PI/3);
            return out;
        }

        static Eigen::Matrix2d oP2D(double scale)
        {
            Eigen::Matrix2d out;
            out << 1,0,
                   0,scale;
            return out;
        }

        static Eigen::Matrix2d oC2D(double scale)
        {
            Eigen::Matrix2d out;
            out << 1,0.5,
                   0,scale;
            return out;
        }
        
        static Eigen::Matrix2d mP2D(double scale, double angle)
        {
            Eigen::Matrix2d out;
            out << 1,scale*std::cos(angle),
                   0,scale*std::sin(angle);
            return out;
        }

        
        static Eigen::Matrix3d cP3D()
        {
            Eigen::Matrix3d out;
            out << 1,0,0,
                   0,1,0,
                   0,0,1;
            return out;
        }

        static Eigen::Matrix3d cI3D()
        {
            Eigen::Matrix3d out;
            out << 1,0,0.5,
                   0,1,0.5,
                   0,0,0.5;
            return out;
        }

        static Eigen::Matrix3d cF3D()
        {
            Eigen::Matrix3d out;
            out << 0,0.5,0.5,
                   0.5,0,0.5,
                   0.5,0.5,0;
            return out;
        }

        static Eigen::Matrix3d hP3D(double height)
        {
            Eigen::Matrix3d out;
            out << 1,std::cos(M_PI/3),0,
                   0,std::sin(M_PI/3),0,
                   0,0,height;
            return out;
        }

        static Eigen::Matrix3d tP3D(double height)
        {
            Eigen::Matrix3d out;
            out << 1,0,0,
                   0,1,0,
                   0,0,height;
            return out;
        }

        static Eigen::Matrix3d tI3D(double height)
        {
            Eigen::Matrix3d out;
            out << 1,0,0.5,
                   0,1,0.5,
                   0,0,height/2;
            return out;
        }

        static Eigen::Matrix3d oP3D(double a,double b, double c)
        {
            Eigen::Matrix3d out;
            out << a,0,0,
                   0,b,0,
                   0,0,c;
            return out;
        }

        static Eigen::Matrix3d oS3D(double a,double b, double c)
        {
            Eigen::Matrix3d out;
            out << a,a/2,0,
                   0,b/2,0,
                   0,0,c;
            return out;
        }
        static Eigen::Matrix3d oI3D(double a,double b, double c)
        {
            Eigen::Matrix3d out;
            out << a,0,a/2,
                   0,b,b/2,
                   0,0,c/2;
            return out;
        }
        
        static Eigen::Matrix3d oF3D(double a,double b, double c)
        {
            Eigen::Matrix3d out;
            out << 0,a/2,a/2,
                   b/2,0,b/2,
                   c/2,c/2,0;
            return out;
        }

        static Eigen::Matrix3d mP3D(double a,double b, double c, double angle)
        {
            Eigen::Matrix3d out;
            out << a,0,c*std::cos(angle),
                   0,b,0,
                   0,0,c*std::sin(angle);
            return out;
        }

        static Eigen::Matrix3d mS3D(double a,double b, double c, double angle)
        {
            Eigen::Matrix3d out;
            out << a,a/2,c*std::cos(angle),
                   0,b/2,0,
                   0,0,c*std::sin(angle);
            return out;
        }

        static Eigen::Matrix3d aP3D(double a,double b, double c, double angle1, double angle2)
        {
            Eigen::Matrix3d out;
            out << a,a*std::cos(angle1),c*std::cos(angle2),
                   0,b*std::cos(angle1),0,
                   0,0,c*std::sin(angle2);
            return out;
        }
        
    };
        //----------------------------------------------------------------------------------------------------------

}
