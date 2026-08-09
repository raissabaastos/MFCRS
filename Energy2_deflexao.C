{
    gROOT->Reset();
    gStyle->SetOptStat(1110);

    FILE *fp = fopen("angulos_deflexao_ROOT.txt", "r");

    if (fp == NULL) {
        Error("Energy2_deflexao", "Nao foi possivel abrir angulos_deflexao_ROOT.txt");
        return;
    }

    const Int_t nEnergies = 5;
    Double_t energyValues[nEnergies] = {
        1.0e16,
        1.0e17,
        1.0e18,
        1.0e19,
        1.0e20
    };

    TH1F *hDeflection[nEnergies];

    for (Int_t i = 0; i < nEnergies; i++) {
        hDeflection[i] = new TH1F(
            Form("hDeflection_%d", i),
            Form("Distribuicao do angulo de deflexao - E = %.0e eV", energyValues[i]),
            90,
            0.0,
            180.0
        );

        hDeflection[i]->GetXaxis()->SetTitle("Angulo de deflexao total [graus]");
        hDeflection[i]->GetYaxis()->SetTitle("Numero de particulas");
    }

    Double_t energy;
    Int_t simulation;
    Double_t deflectionDeg;
    Double_t timeYears;
    Int_t nlines = 0;

    // Le todas as linhas existentes. Nao depende mais de um numero fixo.
    while (fscanf(fp, "%lf %d %lf %lf", &energy, &simulation,
                  &deflectionDeg, &timeYears) == 4) {

        for (Int_t i = 0; i < nEnergies; i++) {
            // Comparacao com tolerancia para valores em ponto flutuante.
            if (TMath::Abs(energy - energyValues[i]) / energyValues[i] < 1.0e-8) {
                hDeflection[i]->Fill(deflectionDeg);
                break;
            }
        }

        nlines++;
    }

    fclose(fp);

    TCanvas *c1 = new TCanvas(
        "c1",
        "Distribuicoes dos angulos de deflexao",
        100,
        100,
        1200,
        800
    );

    c1->Divide(3, 2);

    for (Int_t i = 0; i < nEnergies; i++) {
        c1->cd(i + 1);
        hDeflection[i]->Draw("HIST");
    }

    c1->Update();

    printf("Foram lidas %d particulas do arquivo de deflexoes.\n", nlines);
}
