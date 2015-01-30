import numpy as np
import plotly.plotly as py
from plotly.graph_objs import *

py.sign_in('Python-Demo-Account', 'gwt101uhh0')

polyfit = np.polynomial.polynomial.polyfit
polyval = np.polyval

def qual_ord(phred_string, offset=33):
    return range(len(phred_string)), map(lambda x: (ord(x)-offset), phred_string)

def make_trace(x, y, name=None, mode='markers'):
    trace = Scatter(
        x=x,
        y=y,
        mode=mode, # 'lines+markers'
        #line=Line(shape='spline'),
        name=name
        )
    return trace

def plot_traces(traces, filename='temp_plot3.svg'):
    data = Data(traces)
    fig = Figure(data=data)
    py.image.save_as(fig, filename)

def plot_temp(x, y):
    trace1 = Scatter(x=x, y=y, mode='markers')
    data = Data([trace1])
    layout = Layout(yaxis=YAxis(title='PHRED score'))
    fig = Figure(data=data, layout=layout)
    py.image.save_as(fig, 'temp_plot.svg')

# Don't use anymore.
def fit_and_plot(x, y, degree):
    coeffs = polyfit(x, y, degree)  # watch out! polyval takes these coefficients in reverse order!
    y_fitted = polyval(coeffs[::-1], x)
    # ss = make_trace(x, y_fitted)
    # ss['name'] = 'poly %d' % degree
    orig_trace = Scatter(x=x, y=y, mode='markers')
    fitted_trace = Scatter(x=x, y=y_fitted, mode='markers', name='degree %d' % degree)
    plot_traces([orig_trace, fitted_trace])



#x, y = qual_ord('FEIJJLKLIKJJJFKM?GKJINNHHMKMMMFHLMLHLKMHHJKMLELKLDLJLKL9NLKMLKOKMKJLKEML7IELEIKJJMJLLJIJKIEBIEGDCBA=')

qstrings = '''FGIJJJKFLLGLGPKMFMMILMHMMKIKIKMJKLLKLMHIIDKNHLIKKOLJLNELNLLFKKKKIM8FKLNJJ?KLJIKIKLEILKNJDILJIHDECAA@
FEIFGKFFLLLJGDKKFGHCINNHHIHEJKMHKDLKLH;JLHGMHMJJLKIJLKFLJLCHKKKIIMJGLFNJJIKJJJA?G6HJEMAJK1IB?FGECBA5
BEIGHGKLHHMLGFKMHGKNINNMMHKMJK=JGMIHLHMHKJKILLLINKLJKGELKLKHGKKIBE8GJFJJLLKMJLKJHLE7>KIJK>HG6BDACA3C
EEIBJKKLJLDJILKMKMJKMMNMJMJEJIMHHGJLLKMJLKKLLHNKLOJJBKFLILKMJMKILMJ,KMNLJLKJMIMJKMJHI@LEKM>JFGGDC?AC
F?F6AGCBBBDJBGKM/AHGIGCDMF?@I/=BG;EDHH;CGN=ILHL:LAL;KGFA:LG@-KI7IDJ7K8@F7I7987KJ6GKIMJIJKIIE??,D>+>&
B?IJILCFLGGJJKKMKMMJHKNDMMJKMLIKNGLIHH4ILNKJLMNKLILJ?KKCKLMMKJKIIMJJKKLJ7L7JJJJJDLJJEJNJKM>5,E@AC?A=
FEEBILELJHJJLPCMHIJMHFNMM<EIIIIKKLJHLMMGIKKMLMII:KJJKKK9NLKFH,IIIK8FKMM7JIFLEI6JGIE7'FI76JLBIHD6AA>=
BEGGJKKFLLKFKDKMIGKNANHMHMEIIIMMKLLHHHMHIJKJHHJLLOIJKKL-NLI8:MLJMD8ALKJJ7GKMJIJJJLKJEKIJKMIJJGGACAAA
F9@=5C:JBBLL7:+C:AA:4AN32<2E.<F==.>=L/;:DKKM:L;BLACHDK:--FKFCF'7B7J,6L77J?7C87676G5'5@IJ%>6G?F6,A?3=
AEDBJIKFLKJGGKKG;MKNIMHMMHKMMGIMHHLK=KMLHHLJLENHD:HEMGLGALCMK8LJBEBJKKMIJ?CIEJ6J@>>LJJIJKDHGBEGDA5?=
EGIJIIKLLLJJJMIMKMKIIHHFEMIIM<MMKJGHGMMJLNKLLEKILKCDMLKLHLKMKG7KHKHHKKLFJJKGHJ@JJI>LMK7J6GHEIBGDC<A:
FGDJJKKIKGLKGJKMIMKEHKNMMHJKILFKHMKHLKMNLDKLLHNHKOLHLKKLILMKKKKIBMJHKLNJLICMJJJJNMHHIJLEKI6>IHEDCBA5
EEI=JCCBLKLGJ5KC2AJJFHCDHHKIMCFFKDLH=HFCIKKE:H;I:D:H9:LCAA;>9J,KJ7HJGKEL7JEM8J@J@,HI5FN7KD%7IB,6,??(
HGILJKKJJLKLGKKMKGKMMINMMMKMMLMKKGLLLHMIHKLILMNJLOLJKKKLNLCGJJLKJMJ,KL7@LLKLAGIJK>JJD7IJK5HJ,,=AA<?@
FGGFJIKLLKGLIJKHEGKGL4N?JMJGJKMHNM>KLHMHDKKLFHGBDA:JK:KGHFLK-BLKJDJ,JK'JDGKME,JJ6G+IDJ7JFDIJ6,DA<<A:
EEIGJKCLFLJGJJKMIIAEHHCHJMJIILMKGJLIHKEIHJGLFHLLKKLJLKKLHLGH:MHKIEG'KLJILLFJA>>B@MK7JFF7=16+IFG,<?3A
F.IKJKGLHKLAOMKMEDKCIHNMMMHMJLMBKDLHFKHGHKKIJHJKNGICLNKLILKMLKOKMJJJKFILL?KMJIJJKMJ7DFIJKIEJJ,GECB3(
BCIFJCKBLLGFBGKCEGNJMANFEME@IMHHKMCIHHHJIJFEAMCHKGLJMGLIK9;JHGHFIMILKMIJJGKMMJKJJLJJJMICCG>I66CDCA>-
HEIJJGKJLMJJIMCMFMHMIKNMMMIKMNMMMJGKLKMNLKKMLLKJLKLJLKML:LIMKJOILMJMKKLJHK@MJIIJNMHJMFLEEMIGIHDD6?A=
BEIBHKKIILKKJMKMKIDIMIJMHM?IJNIKGLLILKMFLNKGLHNBLOIJLKKLNLIMKKIFJMJLFKJLLLFCJ@KJKLJJIKLJDMHJJG=DCB;:
A?@BACCIFKDEGDCEFGAEIKCM?I?KDIMBGJJIFHFFHNKIC;CKN.J;B:FADFGBKF7AJAH,?L,7J7E-A7IIH>5HIFAJFD,>'B,6<5A(
BGIJIKKIBHKEBFIJKGKCLKNMJH?KMIMJGMLIFHMGI;KLF;KKLGJDMKLKGLLGKKKKHJJFKMNJFDKGMJKDNI>JJKNJKMHGF?DD<AA=
FEFKJIKLJLKGIMKHKIAIHKCM?MKKE<MMGJIKLD<HIHKIFMJKJKLJLLEL:LLH:FIAMMHGKK7JFIKJJJEINIKAIK7EDM>GI6DAA<A(
=CIFILCJLLKLBJKMKJHJLHNMEMKEMLMMMDLLLKMIHHFELHKHMOIJLKLLNLKHLMHA9KJG>0J7HGCGM@IBJGKIIMNA+@HGI?,DC+A@
FGIJJKKLJKLJLJKMCMJMMMNMJMHIIMMMNHLKLMMNHKKOLMGLJOLJLLLGHLCHKMIIHMJJJKIL7DKMJJJJHLHJJJOJDDH5FGDACA;=
>#IJJLKIJLJGIKKMHMDKHHGMMMIIJNIMMMLKHHMJIKLLLMNLJKJJLKMLKLK@6KLKHMJLKJJLHLKGM@E7KL>IJMLJD5A3?FE=CA?=
HEIKJKKJKLJJJGKMIMJKMMNMMMLKJLMJMMLKLMMIIKNNLMJJNKLJKLLLNLKMJMIILMIJLLLJJKKIJLKJKMJAMJIJKIHHJEGDCB>>
FE@KJGGILLKLJFEM:GKKMNNMMM?MILFHMJLDLHIHHJKOLJJJLILJLKEGJLKJKKIKJMJKGKLLJLKJJJJJNMJLJKIJKMEGBGGEC?;>
<CGFJJKLKLKJFJFMFMKN4KC3?H?@I/HBGC>HFMF:LKG<JLL:LE7AFABIKLK@HGLIRMHAJEEJ7I7GJJK?,>>ED@7,EDEGF6C6@B=0
CEIKJKKLLKJKJMKJHMHMMKNMEMJGMMKMNHGKLMHJLKNMLHLKLOLJKKMLNLLMLJKIREJKLMLLLLKNJIMIGMKEMKL7KMIGJBDECB>(
EEGJHKKILLDJODKMHAJELKNMHMH@MIMHMGCIGMIILNGILMJLL<:JMKKLNLGKKKHFJKJGJLNJLLKIJIJJDLKHJJIJKDLJFBGEC<>0
FGIKJLKLJKKLLMKMHMMJFMNMMMJIJLIMNMILLMHJKNKJLMNJKOLJLKMLNLLMKMKFIMJHKMMJJJKJMLMJNMEEIJLJKAHHIEDEC<AC
;EDFJKKJKMJKILEMKMKMIMNMMMKMJLMMK(GHGHAHIJDI:EIBJ:H;D-LIHF.,(B7'9'8A>#%J7,,JH7,7,,E,D'',#/6+''@D6+A-
FCEKJGKLHHLJJDKMHJKHLMNFMH2GII=BKMCHLDMGIH=IHHGIND:JMGLGJLLM:KKAB,JJ6EJIJKEIJGEJGCJEIKNJD5E+IHDAC5>5
FEI=JIFLFKKJODIMIMMILMDMHKLIJMMMKGIKFHMJIHK<LEIKLKLJ?GKLILKMKFIFMKJL?LMJII7JHGKJK,HEDKLE=MBII?DE<B?-
CGFJJKKBLKIGFPKMKMKCLNC?MMJIJMMKKLLKGHGJHNKOJHC:NKLJFKMLILKHKJIFJMEHLMILF?BIJJKJNLJLJKIEKILJBGEECBA(
<C@GJKKLJLDJOMKMCMDKLHCMMIEKMIMJNMELHHEGLJFMLLIJFILFLGKEJFKHJ8HKLKJKFKELJKKIJJJJ@GJ7LKLJ@IIHIB@ECBA-
EGIGJIKLJLILGLEMKMGKIKNMMMJKMLHHK;>LLHELLJKLLHGLLKLJLKLLKFKMKKIKHKJLGJNLLJKLMLKJKLKHEMNJFGHBFFGDCBA>
FGIJJKKJHLKJJJKMKGGMIKNMMMJIMNIMKJLKEKMJHKKOLLJKLOIJLKKLILLMKBKFLMHLKKNIJLKLJGKINMKEI@NJDJLHFGGDCBA-
EGIFJLKILKDLGJKMHMNJIMNFMMJKMMMJLLLLLMMGINNGLLJLJOJJLKMLNLLKJKOORMJLLMEJLJKCJJJJNMHHJJFEKM>BIGDAAAA@'''

degrees = [2,5,10,20,30]

for i, qstring in enumerate(qstrings.split('\n')):
    x, y = qual_ord(qstring)
    traces = []
    traces.append(make_trace(x, y, name='original'))
    for deg in degrees:
        coeffs = polyfit(x, y, deg)
        y_fitted = polyval(coeffs[::-1], x)
        traces.append(make_trace(x, y_fitted, name='degree %d' % deg, mode='lines'))

    plot_traces(traces, filename=('hiseq2000_%02d.svg' % i))
    break