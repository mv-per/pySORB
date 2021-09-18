

class gplot
{
	public:
		vector<double> xaxis;
		vector<double> yaxis;
		string title = "";
		string xlabel = "";
		string ylabel = "";
		string gnuplot_command_filename = "command.tmp";
		string data_filename = "data.tmp";
		string style = "lt 7 lc 0 w lp";
		string legend = "";
		string grid = "OFF";
		string term = "qt";
		string saveFileName = "";
		string bmargin = "5";
		double xmin = 0;
		double ymin = 0;
		double xmax = 1;
		double ymax= 1;

		void set_xrange(double T_MIN=0, double T_MAX=1){ xmin = T_MIN; xmax = T_MAX; };
		void set_yrange(double T_MIN=0, double T_MAX=1){ ymin = T_MIN; ymax = T_MAX; };
		void set_title(string T=""){ title = T;};
		void set_xlabel(string T=""){ xlabel = T;};
		void set_ylabel(string T=""){ ylabel = T;};
		void set_bmargin(string T="1"){ bmargin = T;};
		void set_legend(string T=""){ legend = T;};
		
		void set_grid(string T="OFF"){ grid = T;};
		void set_term(string T="qt"){ term = T;};
		void set_SavePlot(string F=""){saveFileName = F;};
		


		void set_style( string linetype="lp", string lt="1", string lc = "blue", string lw= "1")
		{ 	
			style = "lt " + lt + " lw " + lw + " lc rgb \"" + lc  + "\" w " + linetype;
		};
		
		void plot(double *x, double *y, size_t *n){
			// size_t n = sizeof(x)/sizeof(double);
			// printf("%d", n);
			std::ofstream fout;
			fout.open(data_filename);
			for (std::size_t i = 0; i < *n; i++){
				fout << x[i] << "\t" << y[i] << endl;
			}
			fout.close();
		};

		void show(){
			std::ofstream gnuplot_command_unit;
			gnuplot_command_unit.open(gnuplot_command_filename.c_str());

			// define margin
			gnuplot_command_unit << "set bmargin "<< bmargin << "\n";
			gnuplot_command_unit << "set title '"<< title << "'\n";
			gnuplot_command_unit << "set key above\n";
			gnuplot_command_unit << "set xlabel '"<< xlabel << "'\n";
    		gnuplot_command_unit << "set ylabel '"<< ylabel << "'\n";
			
			if (xmin){
				gnuplot_command_unit << "set xrange [" << xmin << ":" << xmax << "] \n";
			}
			if (ymin){
				gnuplot_command_unit << "set yrange [" << ymin << ":" << ymax << "] \n";
			}

			if (grid == "ON" || grid == "on"){ 
				gnuplot_command_unit << "set grid\n";
			}
			
			if (term != "qt"){
				gnuplot_command_unit << "set term '"<< term << "'\n";
				if(saveFileName != ""){
					gnuplot_command_unit << "set output " << "\"" << saveFileName << "\"" << "\n";
				}
			}


			gnuplot_command_unit << "plot \"" << data_filename << 
								"\" title '" << legend << "' " << style << " \n";
			gnuplot_command_unit.close();

			// Comando de plot
			string system_command = "gnuplot -persistent " + gnuplot_command_filename;
			system(system_command.c_str());

			//Deleta os arquivos temporários
			remove(gnuplot_command_filename.c_str());
			remove(data_filename.c_str());
		};

		/*    gplot testplot;

    testplot.plot(t,u);
    testplot.set_grid("ON");
    // testplot.set_xrange(1.0,4.0);
    // testplot.set_yrange(0.0,1.0);
    testplot.set_legend("LA124");
    testplot.set_title("LALALALA");
    testplot.set_xlabel("X AXIS");
    testplot.set_ylabel("Y AXIS");
    testplot.set_style("-1", "4", "points");
    // testplot.set_bmargin("5");
    testplot.set_SavePlot("Fig2.png");
    // testplot.set_term("qt");


    testplot.show();
    // std::string xlabel1 = "LALA*/


};