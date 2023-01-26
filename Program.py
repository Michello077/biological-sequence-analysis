import re
import wx
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

'''
Program pobierający lub wczytujący dane w formatach FASTA lub Genbank
i wykonujący na nich podstawowe obliczenia oraz porównania
Program napisali: Michał Humiński i Kacper Kaszuba
studenci 2 roku studiów licencjackich na kierunku Bioinformatyka
w roku akademickim 2020/21
'''

# AB021961.1 - przykładowa sekwencja genu niepodzielonego (ciągłego)
# EU124663.1 - przykładowa sekwencja genu podzielonego (nieciągłego)


def LancuchAminokwasow(evt):
    wybor = lista.GetSelection()
    lancuch = sekwencje[wybor]

    okno1 = wx.Frame(None, title = 'Łańcuch nukleotydów lub aminokwasów', size = (457, 500), pos = (200, 50))
    panel1 = wx.Panel(parent = okno1, pos = (20, 20), size = (400, 400))
    wydruk1 = wx.ListBox(parent = panel1, pos = (20, 20), size = (400, 400))

    while len(lancuch) > 50:
        wydruk1.Append(lancuch[0:49])
        lancuch = lancuch[49:]

    wydruk1.Append(lancuch)
    okno1.Show()


def LiczbaAminokwasow(evt):
    wybor = lista.GetSelection()
    lancuch = sekwencje[wybor]

    aminokwasy = dict()
    for aminokwas in lancuch:
        if aminokwas in aminokwasy.keys():
            aminokwasy[aminokwas] = aminokwasy.get(aminokwas) + 1
        else:
            aminokwasy.update({aminokwas: 1})
    posortowane_aminokwasy = sorted(aminokwasy.items())

    okno2 = wx.Frame(None, title = 'Łańcuch nukleotydów lub aminokwasów', size = (400, 100), pos = (200, 50))
    panel2 = wx.Panel(parent = okno2, pos = (20, 20), size = (200, 100))
    wx.StaticText(panel2, -1, pos = (20, 20), label = f'Liczba poszczególnych nukleotydów lub aminokwasów: \n' + str(posortowane_aminokwasy))

    okno2.Show()


def LiczbaParAminokwasow(evt):
    wybor = lista.GetSelection()
    lancuch = sekwencje[wybor]

    pary_aminokwasow = dict()
    while len(lancuch) > 1:
        para = lancuch[:2]
        if para in pary_aminokwasow.keys():
            pary_aminokwasow[para] = pary_aminokwasow.get(para) + 1
        else:
            pary_aminokwasow.update({para: 1})
        lancuch = lancuch[1:]

    okno3 = wx.Frame(None, title = 'Liczba par aminokwasów', size = (400, 560), pos = (200, 50))
    panel3 = wx.Panel(parent = okno3, pos = (20, 20), size = (400, 400))

    A = 20
    B = 20
    for element in pary_aminokwasow:
        wx.StaticText(panel3, -1, pos = (B, A), label = f'{element}, {pary_aminokwasow[element]}')
        A += 20
        if A == 500:
            B += 50
            A = 20

    okno3.Show()


def LiczbaTrojekNukleotydow(evt):
    wybor = lista.GetSelection()
    lancuch = sekwencje[wybor]

    trojki_nukleotydow = dict()
    while len(lancuch) > 2:
        trojka = lancuch[:3]
        if trojka in trojki_nukleotydow.keys():
            trojki_nukleotydow[trojka] = trojki_nukleotydow.get(trojka) + 1
        else:
            trojki_nukleotydow.update({trojka: 1})
        lancuch = lancuch[1:]

    okno4 = wx.Frame(None, title = 'Liczba trójek nukleotydów', size = (400, 560), pos = (200, 50))
    panel4 = wx.Panel(parent = okno4, pos = (20, 20), size = (400, 400))

    A = 20
    B = 20
    for element in trojki_nukleotydow:
        wx.StaticText(panel4, -1, pos = (B, A), label = f'{element}, {trojki_nukleotydow[element]}')
        A += 20
        if A == 500:
            B += 60
            A = 20

    okno4.Show()


def OtworzFasta(evt):
    dialog = wx.FileDialog(okno, message = 'Wybierz plik FASTA', defaultFile = '', wildcard = '*.FASTA', style = wx.FD_OPEN, pos = (10, 10))

    if dialog.ShowModal() == wx.ID_OK:
        plik = dialog.GetPaths()
        plik_tekst = open(plik[0], 'r')
        odczyt = plik_tekst.read()
        odczyt = odczyt.splitlines()
        odczyt = odczyt[:-1]
        plik_tekst.close()

        sekwencje[:] = []
        nazwy[:] = []
        licznik = 0

        for rekord in odczyt:
            if rekord[0] == '>':
                nazwy.append(rekord)
                sekwencje.append('')
                licznik += 1
            else:
                sekwencje[licznik - 1] += rekord

        lista.InsertItems(nazwy, 0)
        lista.Show()
        dialog.Destroy()


def OtworzGenBank(evt):
    dialog = wx.FileDialog(okno, message = 'Wybierz plik GenBank', defaultFile = '', wildcard = '*.gb', style = wx.FD_OPEN, pos = (10, 10))

    if dialog.ShowModal() == wx.ID_OK:
        plik = dialog.GetPaths()
        plik_tekst = open(plik[0], 'r')
        odczyt = plik_tekst.read().replace('\n', '')

        sekwencje[:] = []
        nazwy[:] = []
        CDS_1[:] = []

        for rekord in SeqIO.parse(plik[0], 'genbank'):
            nazwy.append(rekord.id + ' - sekwencja nukleotydowa')
            nazwy.append(rekord.id + ' - sekwencja kodująca')
            nazwy.append(rekord.id + ' - sekwencja aminokwasowa')

        origin = re.search(r'ORIGIN(.*?)/', odczyt).group(1)
        origin = origin.replace(' ', '')
        origin = ''.join([nukleotyd for nukleotyd in origin if not nukleotyd.isdigit()])
        sekwencje.append(origin)

        CDS = re.findall(r'CDS(.*?)/', odczyt, re.DOTALL)
        laczne_sekwencje_CDS = []
        for lokalizacja_nukleotydow in CDS:
            lokalizacja_nukleotydow_bez_przerw = lokalizacja_nukleotydow.replace(' ', '')
            laczne_sekwencje_CDS.append(lokalizacja_nukleotydow_bez_przerw)
        lokalizacja_sekwencji_CDS = []
        for lokalizacja_CDS in laczne_sekwencje_CDS:
            lokalizacja_CDS_z_numerami = re.findall('[0-9]+', lokalizacja_CDS, re.DOTALL)
            lokalizacja_sekwencji_CDS.extend(lokalizacja_CDS_z_numerami)
        CDS_1.append(lokalizacja_sekwencji_CDS[0])
        sekwencja_kodujaca = ''
        licznik_1 = 0
        licznik_2 = 1
        licznik_3 = 0
        powtarzanie = len(lokalizacja_sekwencji_CDS)
        while licznik_2 < powtarzanie:
            poczatek = int(lokalizacja_sekwencji_CDS[licznik_1]) - 1
            koniec = int(lokalizacja_sekwencji_CDS[licznik_2])
            sekwencja_czastkowa = origin[poczatek:koniec]
            komplementarnosc = laczne_sekwencje_CDS[licznik_3]
            znajdowanie = re.search('complement', komplementarnosc)
            if znajdowanie is None:
                None
            else:
                sekwencja_od_tylu = sekwencja_czastkowa[::-1]
                sekwencja_czastkowa = sekwencja_od_tylu.translate(str.maketrans('tagc', 'atcg'))
            sekwencja_kodujaca += sekwencja_czastkowa
            licznik_1 += 2
            licznik_2 += 2
            licznik_3 += 1
        sekwencje.append(sekwencja_kodujaca)

        translacja = re.findall(r'/translation="(.*?)"', odczyt, re.DOTALL)
        laczna_sekwencja_bialkowa = []
        for sekwencja_aminokwasowa in translacja:
            sekwencja_aminokwasowa_bez_przerw = sekwencja_aminokwasowa.replace(' ', '')
            laczna_sekwencja_bialkowa.append(sekwencja_aminokwasowa_bez_przerw)
        bialko = ''.join(laczna_sekwencja_bialkowa)
        sekwencje.append(bialko)

        plik_tekst.close()
        lista.InsertItems(nazwy, 0)
        lista.Show()
        dialog.Destroy()


def PobierzFasta(evt):
    dialog1 = wx.TextEntryDialog(None, message = 'Proszę podać adres email\nJest potrzebny do pobrania danych z bazy NCBI')
    if dialog1.ShowModal() == wx.ID_OK:
        email = dialog1.GetValue()

        dialog2 = wx.TextEntryDialog(None, message = 'Podaj ID rekordu')
        if dialog2.ShowModal() == wx.ID_OK:
            ID = dialog2.GetValue()

            sekwencje[:] = []
            nazwy[:] = []
            Entrez.email = email
            wyszukanie = Entrez.esearch(db = 'nucleotide', term = ID)
            rekord = Entrez.read(wyszukanie)

            elementy = rekord['IdList']
            for element in elementy:
                wyszukanie = Entrez.efetch(db = 'nucleotide', id = element, rettype = 'FASTA')
                odczyt = wyszukanie.read()
                odczyt = odczyt.splitlines()
                odczyt = odczyt[:-1]
                wyszukanie.close()
                sekwencje[:] = []
                nazwy[:] = []
                licznik = 0
                for rekord in odczyt:
                    if rekord[0] == '>':
                        nazwy.append(rekord)
                        sekwencje.append('')
                        licznik += 1
                    else:
                        sekwencje[licznik - 1] += rekord

                lista.InsertItems(nazwy, 0)
                lista.Show()
            dialog1.Destroy()
            dialog2.Destroy()

        else:
            None

    else:
        None


def PobierzGenBank(evt):
    dialog1 = wx.TextEntryDialog(None, message = 'Proszę podać adres email\nJest potrzebny do pobrania danych z bazy NCBI')
    if dialog1.ShowModal() == wx.ID_OK:
        email = dialog1.GetValue()

        dialog2 = wx.TextEntryDialog(None, message = 'Podaj ID rekordu')
        if dialog2.ShowModal() == wx.ID_OK:
            ID = dialog2.GetValue()

            sekwencje[:] = []
            nazwy[:] = []
            CDS_1[:] = []
            Entrez.email = email
            wyszukanie = Entrez.esearch(db = 'nucleotide', term = ID)
            rekord = Entrez.read(wyszukanie)

            elementy = rekord['IdList']
            for element in elementy:
                wyszukanie = Entrez.efetch(db = 'nucleotide', id = element, rettype = 'gb')
                dane = wyszukanie.read().replace('\n','')

                origin = re.search(r'ORIGIN(.*?)/', dane).group(1)
                origin = origin.replace(' ', '')
                origin = ''.join([nukleotyd for nukleotyd in origin if not nukleotyd.isdigit()])
                nazwy.append(ID + ' - sekwencja nukleotydowa')
                sekwencje.append(origin)

                CDS = re.findall(r'CDS(.*?)/', dane, re.DOTALL)
                laczne_sekwencje_CDS = []
                for lokalizacja_nukleotydow in CDS:
                    lokalizacja_nukleotydow_bez_przerw = lokalizacja_nukleotydow.replace(' ', '')
                    laczne_sekwencje_CDS.append(lokalizacja_nukleotydow_bez_przerw)
                lokalizacja_sekwencji_CDS = []
                for lokalizacja_CDS in laczne_sekwencje_CDS:
                    lokalizacja_CDS_z_numerami = re.findall('[0-9]+', lokalizacja_CDS, re.DOTALL)
                    lokalizacja_sekwencji_CDS.extend(lokalizacja_CDS_z_numerami)
                CDS_1.append(lokalizacja_sekwencji_CDS[0])
                sekwencja_kodujaca = ''
                licznik_1 = 0
                licznik_2 = 1
                licznik_3 = 0
                powtarzanie = len(lokalizacja_sekwencji_CDS)
                while licznik_2 < powtarzanie:
                    poczatek = int(lokalizacja_sekwencji_CDS[licznik_1]) - 1
                    koniec = int(lokalizacja_sekwencji_CDS[licznik_2])
                    sekwencja_czastkowa = origin[poczatek:koniec]
                    komplementarnosc = laczne_sekwencje_CDS[licznik_3]
                    znajdowanie = re.search('complement', komplementarnosc)
                    if znajdowanie is None:
                        None
                    else:
                        sekwencja_od_tylu = sekwencja_czastkowa[::-1]
                        sekwencja_czastkowa = sekwencja_od_tylu.translate(str.maketrans('tagc', 'atcg'))
                    sekwencja_kodujaca += sekwencja_czastkowa
                    licznik_1 += 2
                    licznik_2 += 2
                    licznik_3 += 1
                nazwy.append(ID + ' - sekwencja kodująca')
                sekwencje.append(sekwencja_kodujaca)

                translacja = re.findall(r'/translation="(.*?)"', dane, re.DOTALL)
                laczna_sekwencja_bialkowa = []
                for sekwencja_aminokwasowa in translacja:
                    sekwencja_aminokwasowa_bez_przerw = sekwencja_aminokwasowa.replace(' ', '')
                    laczna_sekwencja_bialkowa.append(sekwencja_aminokwasowa_bez_przerw)
                bialko = ''.join(laczna_sekwencja_bialkowa)
                nazwy.append(ID + ' - sekwencja aminokwasowa')
                sekwencje.append(bialko)

            wyszukanie.close()
            lista.InsertItems(nazwy, 0)
            lista.Show()
            dialog1.Destroy()
            dialog2.Destroy()

        else:
            None

    else:
        None


def PorownaniePrzetlumaczonychSekwencjiOriginiCDS(evt):
    dialog = wx.MessageDialog(None, ('Czy została wybrana sekwencja nukleotydowa?'), ('UWAGA!'), wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION).ShowModal()

    if dialog == wx.ID_YES:
        wybor = lista.GetSelection()
        sekwencja_origin = sekwencje[wybor]
        sekwencja_CDS = Seq(sekwencje[wybor + 1])
        if int(CDS_1[0]) % 3 == 0:
            sekwencja_origin = Seq(sekwencja_origin[2:])
        elif int(CDS_1[0]) % 3 == 1:
            sekwencja_origin = Seq(sekwencja_origin)
        else:
            sekwencja_origin = Seq(sekwencja_origin[1:])
        przetlumaczona_sekwencja_origin = sekwencja_origin.translate()
        przetlumaczona_sekwencja_CDS = sekwencja_CDS.translate()
        dlugosc_origin = len(przetlumaczona_sekwencja_origin)
        dlugosc_CDS = len(przetlumaczona_sekwencja_CDS)

        aminokwasy = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', '*']
        roznice = []
        for aminokwas in aminokwasy:
            roznica = round((przetlumaczona_sekwencja_origin.count(aminokwas)/dlugosc_origin * 100) - (przetlumaczona_sekwencja_CDS.count(aminokwas)/dlugosc_CDS * 100), 2)
            roznice.append(roznica)

        okno5 = wx.Frame(None, title = 'Wynik', size = (300, 600), pos = (200, 50))
        panel5 = wx.Panel(parent = okno5, pos = (20, 20), size = (400, 400))
        wx.StaticText(panel5, pos = (20,20), label = 'Porównanie Przetłumaczonych \nSekwencji ORIGIN i CDS\n')

        szerokosc = 20
        wysokosc = 60
        for element_1 in aminokwasy:
            wx.StaticText(panel5, pos = (szerokosc, wysokosc), label = f'{element_1} :')
            szerokosc += 10
            if szerokosc == 30:
                wysokosc += 20
                szerokosc = 20

        A = 40
        B = 60
        for element_2 in roznice:
            wx.StaticText(panel5, pos = (A, B), label = f'{element_2}%')
            A += 10
            if A == 50:
                B += 20
                A = 40

        okno5.Show()

    else:
        None


def PorownanieSekwencjiCDSibialkowej(evt):
    dialog = wx.MessageDialog(None, ('Czy została wybrana sekwencja kodująca?'), ('UWAGA!'), wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION).ShowModal()

    if dialog == wx.ID_YES:
        wybor = lista.GetSelection()
        sekwencja_CDS = Seq(sekwencje[wybor])
        sekwencja_aminokwasowa = Seq(sekwencje[wybor + 1])
        przetlumaczona_sekwencja = sekwencja_CDS.translate()
        przetlumaczona_sekwencja = str(przetlumaczona_sekwencja)
        przetlumaczona_sekwencja_CDS_lista = []
        for aminokwas in przetlumaczona_sekwencja:
            przetlumaczona_sekwencja_bez_stop = aminokwas.replace('*', '')
            przetlumaczona_sekwencja_CDS_lista.append(przetlumaczona_sekwencja_bez_stop)
        przetlumaczona_sekwencja_CDS = ''.join(przetlumaczona_sekwencja_CDS_lista)

        dlugosc_CDS = len(przetlumaczona_sekwencja_CDS)
        licznik = 0
        licznik_roznic = 0
        roznica_1 = []
        roznica_2 = []
        pozycja = []

        okno6 = wx.Frame(None, title = 'Wynik', size = (400, 600), pos = (200, 50))
        panel6 = wx.Panel(parent = okno6, pos = (20, 20), size = (400, 400))

        if przetlumaczona_sekwencja_CDS == sekwencja_aminokwasowa:
            wx.StaticText(panel6, pos = (20, 20), label = 'Obie sekwencje są identyczne')
        else:
            while dlugosc_CDS > licznik:
                if przetlumaczona_sekwencja_CDS[licznik] == sekwencja_aminokwasowa[licznik]:
                    licznik += 1
                else:
                    roznica_1.append(przetlumaczona_sekwencja_CDS[licznik])
                    roznica_2.append(sekwencja_aminokwasowa[licznik])
                    pozycja.append(licznik)
                    licznik += 1
                    licznik_roznic += 1
            podobienstwo = round(100 - ((licznik_roznic / dlugosc_CDS) * 100), 2)

            wx.StaticText(panel6, pos = (20,20), label = f'Sekwencje są podobne w {podobienstwo}% \nZmiana aminokwasu:')

            S1 = 20
            W1 = 60
            for element_1 in roznica_1:
                wx.StaticText(panel6, pos = (S1,W1), label = f'{element_1}')
                S1 += 10
                if S1 == 30:
                    W1 += 20
                    S1 = 20

            S2 = 40
            W2 = 60
            for element_2 in roznica_2:
                wx.StaticText(panel6, pos = (S2, W2), label = f'na {element_2}')
                S2 += 10
                if S2 == 50:
                    W2 += 20
                    S2 = 40

            S3 = 80
            W3 = 60
            for element_3 in pozycja:
                wx.StaticText(panel6, pos = (S3, W3), label = f'w pozycji {element_3}')
                S3 += 15
                if S3 == 95:
                    W3 += 20
                    S3 = 80

        okno6.Show()

    else:
        None


def PorownanieSekwencjiOriginiCDS(evt):
    dialog = wx.MessageDialog(None, ('Czy została wybrana sekwencja nukleotydowa?'), ('UWAGA!'), wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION).ShowModal()

    if dialog == wx.ID_YES:
        wybor = lista.GetSelection()
        sekwencja_origin = sekwencje[wybor]
        sekwencja_CDS = sekwencje[wybor + 1]
        dlugosc_origin = len(sekwencja_origin)
        dlugosc_CDS = len(sekwencja_CDS)

        roznica_a = round((sekwencja_origin.count('a')/dlugosc_origin * 100) - (sekwencja_CDS.count('a')/dlugosc_CDS * 100), 2)
        roznica_t = round((sekwencja_origin.count('t')/dlugosc_origin * 100) - (sekwencja_CDS.count('t')/dlugosc_CDS * 100), 2)
        roznica_c = round((sekwencja_origin.count('c')/dlugosc_origin * 100) - (sekwencja_CDS.count('c')/dlugosc_CDS * 100), 2)
        roznica_g = round((sekwencja_origin.count('g')/dlugosc_origin * 100) - (sekwencja_CDS.count('g')/dlugosc_CDS * 100), 2)

        okno7 = wx.Frame(None, title = 'Wynik', size = (400, 200), pos = (200, 50))
        panel7 = wx.Panel(parent = okno7, pos = (20, 20), size = (400, 400))
        wx.StaticText(panel7, pos = (20, 20),
                      label=f'Procentowa różnica zawartości nukleotydów \npomiędzy sekwencjami ORIGIN i CDS: \n'
                            f'\nAdenina: {roznica_a}%\nTymina: {roznica_t}%\nCytozyna: {roznica_c}%\nGuanina: {roznica_g}%')
        okno7.Show()

    else:
        None


def Zamknij(evt):
    dialog = wx.MessageDialog(okno, 'Czy na pewno?', 'Kończymy pracę', style = wx.OK | wx.CANCEL)
    wyjscie = dialog.Showmodal()
    dialog.Destroy()
    if wyjscie == wx.ID_OK:
        okno.Close()


Program = wx.App()
sekwencje = []
nazwy = []
CDS_1 = []

okno = wx.Frame(None, title = 'Menu programu', size = (800, 600), pos = (200, 50))

menulistwa = wx.MenuBar()

menu1 = wx.Menu()
menu1_opcja1 = menu1.Append(wx.ID_ANY, 'Wczytaj Fasta', 'Czytaj Dane')
okno.Bind(wx.EVT_MENU, OtworzFasta, menu1_opcja1)
menu1_opcja2 = menu1.Append(wx.ID_ANY, 'Wczytaj GenBank', 'Czytaj Dane')
okno.Bind(wx.EVT_MENU, OtworzGenBank, menu1_opcja2)
menu1_opcja3 = menu1.Append(wx.ID_ANY, 'Pobierz Fasta', 'Czytaj Dane')
okno.Bind(wx.EVT_MENU, PobierzFasta, menu1_opcja3)
menu1_opcja4 = menu1.Append(wx.ID_ANY, 'Pobierz GenBank', 'Czytaj Dane')
okno.Bind(wx.EVT_MENU, PobierzGenBank, menu1_opcja4)
menulistwa.Append(menu1, 'Dane')

menu2 = wx.Menu()
menu2_opcja1 = menu2.Append(wx.ID_ANY, 'Łańcuch nukleotydów/aminokwasów', 'Oblicz')
okno.Bind(wx.EVT_MENU, LancuchAminokwasow, menu2_opcja1)
menu2_opcja2 = menu2.Append(wx.ID_ANY, 'Liczba nukleotydów/aminokwasów', 'Oblicz')
okno.Bind(wx.EVT_MENU, LiczbaAminokwasow, menu2_opcja2)
menu2_opcja3 = menu2.Append(wx.ID_ANY, 'Liczba par aminokwasów', 'Oblicz')
okno.Bind(wx.EVT_MENU, LiczbaParAminokwasow, menu2_opcja3)
menu2_opcja4 = menu2.Append(wx.ID_ANY, 'Liczba trójek nukleotydów', 'Oblicz')
okno.Bind(wx.EVT_MENU, LiczbaTrojekNukleotydow, menu2_opcja4)
menulistwa.Append(menu2, 'Obliczenia')

menu3 = wx.Menu()
menu3_opcja1 = menu3.Append(wx.ID_ANY, 'Porównanie przetłumaczonej sekwencji CDS i aminokwasowej' ,'Porównaj')
okno.Bind(wx.EVT_MENU, PorownanieSekwencjiCDSibialkowej, menu3_opcja1)
menu3_opcja2 = menu3.Append(wx.ID_ANY, 'Porównanie sekwencji nukleotydowych origin i CDS' ,'Porównaj')
okno.Bind(wx.EVT_MENU, PorownanieSekwencjiOriginiCDS, menu3_opcja2)
menu3_opcja3 = menu3.Append(wx.ID_ANY, 'Porównanie przetłumaczonych sekwencji aminokwasowych origin i CDS' ,'Porównaj')
okno.Bind(wx.EVT_MENU, PorownaniePrzetlumaczonychSekwencjiOriginiCDS, menu3_opcja3)
menulistwa.Append(menu3, 'Porównanie sekwencji GenBank')

menu4 = wx.Menu()
menu4_opcja1 = menu4.Append(wx.ID_ANY, 'Koniec', 'Koniec programu')
okno.Bind(wx.EVT_MENU, Zamknij, menu4_opcja1)
menulistwa.Append(menu4, 'Wyjście')

okno.SetMenuBar(menulistwa)

panel = wx.Panel(parent = okno)
lista = wx.ListBox(parent = panel, pos = (20, 20), size = (400, 400))
lista.Hide()

okno.Show()

Program.MainLoop()
