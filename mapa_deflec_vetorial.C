#include <TCanvas.h>
#include <TPolyLine.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TAxis.h>
#include <TH2F.h>
#include <TLine.h>
#include <TLatex.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

struct Ponto {
   double lat;        // radianos
   double lon;        // radianos
   int tipo;          // 0 = escape, 1 = inicial
   double energy;
   int sim;
   double time_years;
};

double normLonRad(double lon)
{
   while (lon > TMath::Pi()) lon -= 2.0 * TMath::Pi();
   while (lon < -TMath::Pi()) lon += 2.0 * TMath::Pi();
   return lon;
}

// Projecao Aitoff manual.
// Entrada: lon, lat em radianos.
// Saida: x, y no plano da projecao.
void Aitoff(double lon, double lat, double &x, double &y)
{
   lon = normLonRad(lon);

   double alpha = lon / 2.0;
   double c = cos(lat) * cos(alpha);
   double denom = sqrt(1.0 + c);

   if (denom < 1e-12) {
      x = 0.0;
      y = 0.0;
      return;
   }

   x = 2.0 * sqrt(2.0) * cos(lat) * sin(alpha) / denom;
   y = sqrt(2.0) * sin(lat) / denom;
}

// Conversao esferica para cartesiano
void SphToCart(double lat, double lon, double &x, double &y, double &z)
{
   x = cos(lat) * cos(lon);
   y = cos(lat) * sin(lon);
   z = sin(lat);
}

// Conversao cartesiano para esferico
void CartToSph(double x, double y, double z, double &lat, double &lon)
{
   double r = sqrt(x*x + y*y + z*z);

   x /= r;
   y /= r;
   z /= r;

   lat = asin(z);
   lon = atan2(y, x);
}

// Desenha uma curva na esfera entre dois pontos.
// Isso evita a reta artificial no plano.
void DrawGreatCirclePath(
   double lat1,
   double lon1,
   double lat2,
   double lon2,
   int color,
   int nsteps = 40
)
{
   double x1, y1, z1;
   double x2, y2, z2;

   SphToCart(lat1, lon1, x1, y1, z1);
   SphToCart(lat2, lon2, x2, y2, z2);

   double dot = x1*x2 + y1*y2 + z1*z2;

   if (dot > 1.0) dot = 1.0;
   if (dot < -1.0) dot = -1.0;

   double omega = acos(dot);

   std::vector<double> xs;
   std::vector<double> ys;

   for (int i = 0; i <= nsteps; i++) {
      double t = double(i) / double(nsteps);

      double x, y, z;

      if (fabs(omega) < 1e-12) {
         x = x1;
         y = y1;
         z = z1;
      } else {
         double a = sin((1.0 - t) * omega) / sin(omega);
         double b = sin(t * omega) / sin(omega);

         x = a*x1 + b*x2;
         y = a*y1 + b*y2;
         z = a*z1 + b*z2;
      }

      double lat, lon;
      CartToSph(x, y, z, lat, lon);

      double xp, yp;
      Aitoff(lon, lat, xp, yp);

      xs.push_back(xp);
      ys.push_back(yp);
   }

   TPolyLine *pl = new TPolyLine(xs.size(), xs.data(), ys.data());
   pl->SetLineColor(color);
   pl->SetLineWidth(1);
   pl->Draw("L");
}

// Desenha a borda da projecao Aitoff
void DrawAitoffBoundary()
{
   const int n = 400;
   double xs[n + 1];
   double ys[n + 1];

   for (int i = 0; i <= n; i++) {
      double t = 2.0 * TMath::Pi() * double(i) / double(n);

      xs[i] = 2.0 * sqrt(2.0) * cos(t);
      ys[i] = sqrt(2.0) * sin(t);
   }

   TPolyLine *border = new TPolyLine(n + 1, xs, ys);
   border->SetLineColor(kBlack);
   border->SetLineWidth(2);
   border->Draw("L");
}

// Desenha grade de latitude e longitude
void DrawGrid()
{
   // Longitudes
   for (int lon_deg = -150; lon_deg <= 150; lon_deg += 30) {
      std::vector<double> xs;
      std::vector<double> ys;

      double lon = lon_deg * TMath::Pi() / 180.0;

      for (int lat_deg = -90; lat_deg <= 90; lat_deg++) {
         double lat = lat_deg * TMath::Pi() / 180.0;

         double x, y;
         Aitoff(lon, lat, x, y);

         xs.push_back(x);
         ys.push_back(y);
      }

      TPolyLine *line = new TPolyLine(xs.size(), xs.data(), ys.data());
      line->SetLineColor(kGray + 1);
      line->SetLineStyle(3);
      line->Draw("L");
   }

   // Latitudes
   for (int lat_deg = -60; lat_deg <= 60; lat_deg += 30) {
      std::vector<double> xs;
      std::vector<double> ys;

      double lat = lat_deg * TMath::Pi() / 180.0;

      for (int lon_deg = -180; lon_deg <= 180; lon_deg++) {
         double lon = lon_deg * TMath::Pi() / 180.0;

         double x, y;
         Aitoff(lon, lat, x, y);

         xs.push_back(x);
         ys.push_back(y);
      }

      TPolyLine *line = new TPolyLine(xs.size(), xs.data(), ys.data());
      line->SetLineColor(kGray + 1);
      line->SetLineStyle(3);
      line->Draw("L");
   }
}

void DrawReferenceLabels()
{
    TLatex *label = new TLatex();

    label->SetTextSize(0.035);
    label->SetTextAlign(22);

    // Longitude
    label->DrawLatex(-2.72, -0.10, "-180^{#circ}");
    label->DrawLatex( 2.72, -0.10, "180^{#circ}");

    // Latitude
    label->DrawLatex(0.0,  1.47, "90^{#circ}");
    label->DrawLatex(0.0, -1.47, "-90^{#circ}");
}


TCanvas *c1 = new TCanvas(
    "c1",
    "Mapas vetoriais de deflexao",
    1200,
    850
);

c1->Divide(2, 2);

   std::ifstream in("lat_long_velocidades_ROOT_ASS_P.csv");

   if (!in.is_open()) {
      std::cerr << "Erro: nao consegui abrir lat_long_velocidades_ROOT_2.csv" << std::endl;
      return c1;
   }

   std::string line;
   std::getline(in, line); // pula o cabecalho

   std::vector<Ponto> iniciais;
   std::vector<Ponto> escapes;

   while (std::getline(in, line)) {
      if (line.empty()) continue;

      std::stringstream ss(line);
      std::string item;

      Ponto p;

      std::getline(ss, item, ',');
      p.lat = std::stod(item);

      std::getline(ss, item, ',');
      p.lon = std::stod(item);

      std::getline(ss, item, ',');
      p.tipo = std::stoi(item);

      std::getline(ss, item, ',');
      p.energy = std::stod(item);

      std::getline(ss, item, ',');
      p.sim = std::stoi(item);

      std::getline(ss, item, ',');
      p.time_years = std::stod(item);

      if (p.tipo == 1) iniciais.push_back(p);
      if (p.tipo == 0) escapes.push_back(p);
   }

   in.close();

   std::cout << "Iniciais: " << iniciais.size() << std::endl;
   std::cout << "Escapes: " << escapes.size() << std::endl;

double energies[4] = {
    1e17,
    1e18,
    1e19,
    1e20
};

int exponents[4] = {
    17,
    18,
    19,
    20
};


for (int i = 0; i < 4; i++) {

    c1->cd(i + 1);

    double targetEnergy = energies[i];


    // Cada histograma precisa ter um nome diferente no ROOT
    TH2F *frame = new TH2F(
        Form("frame_%d", i),
        Form(
            "Mapa de deflex#tilde{a}o (ASS) - pr#acute{o}ton - E = 10^{%d} eV",
            exponents[i]
        ),
        10, -3.0, 3.0,
        10, -1.6, 1.6
    );


    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetLabelSize(0);

    frame->GetXaxis()->SetTickLength(0);
    frame->GetYaxis()->SetTickLength(0);

    frame->Draw();


    DrawGrid();
    DrawAitoffBoundary();
    DrawReferenceLabels();


    // Agora desenha APENAS partículas desta energia
    for (auto &ini : iniciais) {

        bool energia_desejada =
            fabs(ini.energy - targetEnergy) / targetEnergy < 1e-9;

        if (!energia_desejada)
            continue;


        for (auto &esc : escapes) {

            bool mesmo_sim =
                ini.sim == esc.sim;

            bool mesma_energy =
                fabs(ini.energy - esc.energy) / ini.energy < 1e-9;

            bool mesmo_tempo =
                fabs(ini.time_years - esc.time_years) < 1e-6;


            if (
                mesmo_sim &&
                mesma_energy &&
                mesmo_tempo
            ) {

                DrawGreatCirclePath(
                    ini.lat,
                    ini.lon,
                    esc.lat,
                    esc.lon,
                    kMagenta + 2,
                    5
                );

                break;
            }
        }
    }
}

   c1->Update();

   return c1;
}
