using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Text.RegularExpressions;
using System.Diagnostics;

namespace Lenarda_Jonesa
{
    class Data
    {
        public string name { get; set; }
        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }
    }

    class VelocityVector
    {
        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }
    }

    class LocationVector
    {
        public string name { get; set; }
        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }
    }

    class ForceVector
    {
        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }
    }

    class AccelerationVector
    {
        public double x { get; set; }
        public double y { get; set; }
        public double z { get; set; }
    }

    class stage2Data
    {
        public LocationVector _locationVector { get; set; }
        public VelocityVector _velocityVector { get; set; }
        public ForceVector _forceVector { get; set; }
        public AccelerationVector _accelerationVectorT { get; set; }
        public AccelerationVector _accelerationVectorT_H { get; set; }
    }

    class Program
    {
        //data definitions
        //extensions - stage 1
        //string _inputFile = args[0];
        public static string[] columns = new string[12];
        static string _configFileName = string.Empty;
        static double _mass = 0.0;
        //potential
        static double _eps = 0.0;
        static double _sigma = 0.0;
        static double _timestep = 0.0;
        static int _steps = 0;
        static string _thermoFilename = string.Empty;
        static string _trajFilename = string.Empty;

        //zad5 - new variables
        static double _nThermo = 0;
        static double _nTraj = 0;
        static double _rCutoff = 0.0;
        //-------------

        
        static double _mInverted = 0.0;
        static double A = 0.0;
        static double _AInverted = 0.0;


        static int iterator = 0;

        static Data inputDataComponent = new Data();
        static Data Rij = new Data();
        static double ordinaryF = new double();

        static double RijCube = 0;
        static double RijCubeDivided = 0;
        static double RijCubeDividedAndMultipliedBySigmaCube = 0;
        static double RijDivided6 = 0;
        static double RijDivided12 = 0;

        static double factorF = 0.0;
        static double factorV = 0.0;
        static double sigmaToSecondPower = 0.0;
        static double _timestep2 = 0.0;
        //stage2
        static int atomCount;// = 0;

        //count of atoms
        static int lineCount;// = File.ReadLines(_configFileName).Count();
        static stage2Data[] data;// = new stage2Data[lineCount];

        //static Data[] inputDataTable;//= new Data[lineCount];

        static double[] V;// 
        static double _summaryEnergy = 0.0;
        static double _summaryV = 0.0;
        static double _summaryK = 0.0;
        static double _time = 0.0;

        //timers
        static Stopwatch _totalTime = new Stopwatch();
        static Stopwatch _stepTime = new Stopwatch();
        static TimeSpan _totalTimeSpan;
        static TimeSpan _totalStepTimeSpan;

        static void readInputFile(string _inputFile)
        {
            try
            {
                if (!File.Exists(_inputFile))
                    Console.WriteLine("File not existst!");
                else
                {
                    StreamReader strRead = new StreamReader(_inputFile);
                    while (!strRead.EndOfStream)
                    {

                        columns = strRead.ReadLine().Split((char[])null, StringSplitOptions.RemoveEmptyEntries);

                        if ("config".Equals(columns[0]))
                            _configFileName = columns[1];

                        if ("mass".Equals(columns[0]))
                            _mass = Math.Abs(Convert.ToDouble(columns[1]));

                        if ("potential".Equals(columns[0]))
                        {
                            _eps = Math.Abs(Convert.ToDouble(columns[1]));
                            _sigma = Math.Abs(Convert.ToDouble(columns[2]));
                        }

                        if ("timestep".Equals(columns[0]))
                            _timestep = Math.Abs(Convert.ToDouble(columns[1]));

                        if ("steps".Equals(columns[0]))
                            _steps = Math.Abs(Convert.ToInt32(columns[1]));

                        if ("thermo".Equals(columns[0]))
                            _thermoFilename = columns[1];

                        if ("print".Equals(columns[0]))
                        {
                            _nThermo = Math.Abs(Convert.ToDouble(columns[1]));
                            _nTraj = Math.Abs(Convert.ToDouble(columns[2]));
                        }

                        if ("traj".Equals(columns[0]))
                            _trajFilename = columns[1];

                        if("r_cutoff".Equals(columns[0]))
                            _rCutoff = Math.Abs(Convert.ToDouble(columns[1]));
                    }
                    strRead.Close();

                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }
        }

        static void printRunParameters()
        {
            Console.WriteLine(_configFileName);
            Console.WriteLine(_mass);
            Console.WriteLine(_eps);
            Console.WriteLine(_sigma);
            Console.WriteLine(_timestep);
            Console.WriteLine(_steps);
            Console.WriteLine(_thermoFilename);

            Console.WriteLine(_nThermo);
            Console.WriteLine(_nTraj);

            Console.WriteLine(_trajFilename);
        }

        static void readInitialConfiguration()
        {
            try
            {
                lineCount = File.ReadLines(_configFileName).Count();
                data = new stage2Data[lineCount];
                //memory allocation - stage 2
                for (int k = 0; k < lineCount; k++)
                {

                    data[k] = new stage2Data();
                    data[k]._accelerationVectorT = new AccelerationVector();
                    //set acceleration vector of all atoms
                    data[k]._accelerationVectorT.x = 0.0;
                    data[k]._accelerationVectorT.y = 0.0;
                    data[k]._accelerationVectorT.z = 0.0;
                    data[k]._accelerationVectorT_H = new AccelerationVector();
                    //set acceleration vector of all atoms
                    data[k]._accelerationVectorT_H.x = 0.0;
                    data[k]._accelerationVectorT_H.y = 0.0;
                    data[k]._accelerationVectorT_H.z = 0.0;
                    data[k]._forceVector = new ForceVector();
                    data[k]._locationVector = new LocationVector();
                    data[k]._velocityVector = new VelocityVector();
                    //set start speed of all atoms
                    data[k]._velocityVector.x = 0.0;
                    data[k]._velocityVector.y = 0.0;
                    data[k]._velocityVector.z = 0.0;
                }

                if (!File.Exists(_configFileName))
                    Console.WriteLine("File not existst!");
                else
                {
                    StreamReader strRead = new StreamReader(_configFileName);
                    columns = strRead.ReadLine().Split(' ');
                    atomCount = Convert.ToInt32(columns[0]);
                    columns = strRead.ReadLine().Split(' ');
                    while (!strRead.EndOfStream)
                    {
                        columns = strRead.ReadLine().Split((char[])null, StringSplitOptions.RemoveEmptyEntries);

                        //inputDataTable[iterator].name = columns[0];
                        //inputDataTable[iterator].x = Convert.ToDouble(columns[1]);
                        //inputDataTable[iterator].y = Convert.ToDouble(columns[2]);
                        //inputDataTable[iterator].z = Convert.ToDouble(columns[3]);
                        data[iterator]._locationVector.name = columns[0];
                        data[iterator]._locationVector.x = Convert.ToDouble(columns[1]);
                        data[iterator]._locationVector.y = Convert.ToDouble(columns[2]);
                        data[iterator]._locationVector.z = Convert.ToDouble(columns[3]);

                        iterator++;
                    }
                    strRead.Close();
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }
        }

        static void initializeComputations()
        {
            try
            {
                factorF = 24d * _eps;
                factorV = 2d * _eps;
                sigmaToSecondPower = _sigma * _sigma;

                _timestep2 = _timestep * _timestep;
                _mInverted = 1.0 / _mass;
                A = 0.009648545946;
                _AInverted = 1.0 / A;
                //z uwagi na zasade otwierania strumienia tylko w chwili zapisu, 
                //nie otwieram tutaj plikow wyjsciowych
                //tylko robie to tuz przed zapisem danych, 
                //tym sposobem nie mam otwartego niepotrzebnie strumienia zapisu
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }


        }

        static void calculateForces()
        {
            try
            {
                V = new double[atomCount]; //potencjal
                for (int i = 0; i < atomCount; i++)
                {
                    //F[i] = new Data();
                    V[i] = 0.0;
                    data[i]._forceVector.x = 0.0;
                    data[i]._forceVector.y = 0.0;
                    data[i]._forceVector.z = 0.0;
                }

                for (int i = 0; i < atomCount; i++)
                {
                    for (int j = i + 1; j < atomCount; j++)
                    {
                        //Rij - delta x/y/z
                        Rij.x = data[j]._locationVector.x - data[i]._locationVector.x;
                        Rij.y = data[j]._locationVector.y - data[i]._locationVector.y;
                        Rij.z = data[j]._locationVector.z - data[i]._locationVector.z;

                        if (Math.Sqrt(Rij.x * Rij.x + Rij.y * Rij.y + Rij.z * Rij.z) <= _rCutoff)
                        {


                            //do kwadratu
                            RijCube = Rij.x * Rij.x + Rij.y * Rij.y + Rij.z * Rij.z;
                            RijCubeDivided = 1 / RijCube;
                            RijCubeDividedAndMultipliedBySigmaCube = sigmaToSecondPower * RijCubeDivided;

                            RijDivided6 = RijCubeDividedAndMultipliedBySigmaCube * RijCubeDividedAndMultipliedBySigmaCube * RijCubeDividedAndMultipliedBySigmaCube;
                            RijDivided12 = RijDivided6 * RijDivided6;
                            ordinaryF = factorF * RijCubeDivided * (2 * RijDivided12 - RijDivided6);

                            //potencjaly i-ty / j-ty
                            V[i] += (RijDivided12 - RijDivided6);
                            V[j] += (RijDivided12 - RijDivided6);
                        }
                        else
                        {
                            ordinaryF = 0;
                            V[i] += 0;
                            V[j] += 0;
                        }

                        data[i]._forceVector.x -= ordinaryF * Rij.x;
                        data[i]._forceVector.y -= ordinaryF * Rij.y;
                        data[i]._forceVector.z -= ordinaryF * Rij.z;

                        data[j]._forceVector.x += ordinaryF * Rij.x;
                        data[j]._forceVector.y += ordinaryF * Rij.y;
                        data[j]._forceVector.z += ordinaryF * Rij.z;
                    }
                }
                for (int i = 0; i < atomCount; i++)
                {
                    //data[i]._forceVector.x *= factorF;
                    //data[i]._forceVector.y *= factorF;
                    //data[i]._forceVector.z *= factorF;
                    V[i] *= factorV;
                }

            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }
        }

        static void calculateAccelerations()
        {
            try
            {
                for (int i = 0; i < atomCount; i++)
                {
                    data[i]._accelerationVectorT.x = _mInverted * data[i]._forceVector.x * A;
                    data[i]._accelerationVectorT.y = _mInverted * data[i]._forceVector.y * A;
                    data[i]._accelerationVectorT.z = _mInverted * data[i]._forceVector.z * A;
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }

        }

        static void updatePositions()
        {
            try
            {

                for (int i = 0; i < atomCount; i++)
                {
                    data[i]._locationVector.x +=
                        data[i]._velocityVector.x
                        * _timestep + 0.5
                        * data[i]._accelerationVectorT.x * _timestep2;

                    data[i]._locationVector.y +=
                        data[i]._velocityVector.y
                        * _timestep + 0.5
                        * data[i]._accelerationVectorT.y * _timestep2;

                    data[i]._locationVector.z +=
                        data[i]._velocityVector.z
                        * _timestep + 0.5
                        * data[i]._accelerationVectorT.z * _timestep2;
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }

        }

        static void updateAccelerations()
        {
            try
            {
                for (int i = 0; i < atomCount; i++)
                {
                    data[i]._accelerationVectorT_H.x = data[i]._accelerationVectorT.x;
                    data[i]._accelerationVectorT_H.y = data[i]._accelerationVectorT.y;
                    data[i]._accelerationVectorT_H.z = data[i]._accelerationVectorT.z;
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }

        }

        static void updateVelocities()
        {
            try
            {
                for (int i = 0; i < atomCount; i++)
                {
                    //data[i]._velocityVector.x = 0.0;
                    //data[i]._velocityVector.y = 0.0;
                    //data[i]._velocityVector.z = 0.0;
                }

                for (int i = 0; i < atomCount; i++)
                {
                    data[i]._velocityVector.x +=
                        0.5 * _timestep
                        * (data[i]._accelerationVectorT_H.x + data[i]._accelerationVectorT.x);

                    data[i]._velocityVector.y +=
                        0.5 * _timestep
                        * (data[i]._accelerationVectorT_H.y + data[i]._accelerationVectorT.y);

                    data[i]._velocityVector.z +=
                        0.5 * _timestep
                        * (data[i]._accelerationVectorT_H.z + data[i]._accelerationVectorT.z);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }
        }

        static void computeThermo()
        {
            try
            {
                _summaryV = 0.0;
                _summaryK = 0.0;
                _summaryEnergy = 0.0;
                for (int i = 0; i < atomCount; i++)
                {
                    _summaryV += V[i];
                    _summaryK += 0.5 * _mass * (data[i]._velocityVector.x * data[i]._velocityVector.x
                        + data[i]._velocityVector.y * data[i]._velocityVector.y
                        + data[i]._velocityVector.z * data[i]._velocityVector.z
                        );

                }
                _summaryK *= _AInverted;
                _summaryEnergy = _summaryV + _summaryK;
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }
        }



        static void Main(string[] args)
        {
            System.Globalization.CultureInfo customCulture = (System.Globalization.CultureInfo)System.Threading.Thread.CurrentThread.CurrentCulture.Clone();
            customCulture.NumberFormat.NumberDecimalSeparator = ".";
            System.Threading.Thread.CurrentThread.CurrentCulture = customCulture;

            try
            {
                if ("1".Equals(args.Length.ToString()))
                {



                    readInputFile(args[0]);
                    printRunParameters();

                    readInitialConfiguration();
                    initializeComputations();

                    //here timers starts
                    _totalTime.Start();

                    calculateForces();
                    calculateAccelerations();

                    if (!File.Exists(_thermoFilename))
                        File.Create(_thermoFilename);

                    if (!File.Exists(_trajFilename) && !string.Empty.Equals(_trajFilename))
                        File.Create(_trajFilename);

                    if (!string.Empty.Equals(_thermoFilename) && !string.Empty.Equals(_trajFilename))
                    {


                        System.IO.File.WriteAllText(_thermoFilename, string.Empty);
                        System.IO.File.WriteAllText(_trajFilename, string.Empty);

                        System.IO.StreamWriter _thermo = new System.IO.StreamWriter(_thermoFilename, true);
                        System.IO.StreamWriter _thraj = new System.IO.StreamWriter(_trajFilename, true);

                        _thermo.WriteLine(String.Format("{0,20}{1,20}{2,20}{3,20}{4,20}",
                            "# 1 - step",
                            " 2 - time",
                            " 3 - E_tot",
                            " 4 - V_tot",
                            " 5 - K_tot"));

                        Console.WriteLine("Starting..");
                        for (int i = 0; i < _steps; i++)
                        {
                            _time = i * _timestep;

                            _stepTime.Reset();
                            _stepTime.Start();

                            updatePositions();
                            updateAccelerations();

                            calculateForces();
                            calculateAccelerations();
                            updateVelocities();

                            computeThermo();

                            _stepTime.Stop();
                            _totalStepTimeSpan = _stepTime.Elapsed;

                            if (_nTraj != 0 || _nThermo != 0)
                            {
                                if (i % _nThermo == 0)
                                {
                                    //printThermoOnScreen
                                    Console.WriteLine(
                                        i + " "
                                        + _time
                                        + " "
                                        + _summaryEnergy
                                        + " "
                                        + _summaryV + " "
                                        + _summaryK);

                                    //writeThermoToFile - poprawic
                                    _thermo.WriteLine(String.Format("{0,20}{1,20}{2,25}{3,25}{4,25}",
                                    i,
                                    _time,
                                    _summaryEnergy,
                                    _summaryV,
                                    _summaryK
                                    ));
                                }

                                if (i % _nTraj == 0)
                                {
                                    //writeTraj to file
                                    _thraj.WriteLine(atomCount);
                                    _thraj.WriteLine("timestep " + i + ", time " + _time);
                                    for (int j = 0; j < atomCount; j++)
                                    {
                                        _thraj.WriteLine(data[j]._locationVector.name + " " +
                                            data[j]._locationVector.x + " " +
                                            data[j]._locationVector.y + " " +
                                            data[j]._locationVector.z
                                            );
                                    }
                                }
                                
                            }


                            Console.WriteLine("# total time (per step, in ms): " + _totalStepTimeSpan.TotalMilliseconds);
                        }
                        _totalTime.Stop();
                        _thermo.Close();
                        _thraj.Close();
                        _totalTimeSpan = _totalTime.Elapsed;
                        Console.WriteLine("Calculations complete.");
                        Console.WriteLine("# number of steps: " + _steps);
                        Console.WriteLine("# total time (in ms): " + _totalTimeSpan.TotalMilliseconds);
                        

                    }
                    
                    Console.ReadKey();
                }
                else
                {
                    System.Environment.Exit(0);
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.ToString());
            }
        }
    }
}
