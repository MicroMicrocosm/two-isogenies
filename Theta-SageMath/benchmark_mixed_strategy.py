from sage.all import GF, EllipticCurve
from isogeny_diamond import generate_splitting_kernel, DIAMONDS
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny
from theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt
from utilities.supersingular import torsion_basis
from utilities.strategy import optimised_strategy
from utilities.utils import speed_up_sagemath


import time

speed_up_sagemath()


def time_theta(test_index, test_N=1):
    """
    Selects a kernel generating an isogeny between elliptic products from
    `isogeny_diamond.py` using `generate_splitting_kernel()` and times
    the average computation time of codomain and evaluation time for the
    isogeny chain in the theta model.
    """
    (P1, Q1, P2, Q2), _ = generate_splitting_kernel(test_index)

    # Create kernel from CouplePoint data
    ker_Phi = (CouplePoint(P1, P2), CouplePoint(Q1, Q2))
    EA, EB = ker_Phi[0].curves()

    _, ea, eb, _, _ = DIAMONDS[test_index]
    PB3, _ = torsion_basis(EB, 3**eb)

    # Optimal strategy for the isogeny
    strategy = STRATEGY[ea]

    # start profiler
    t0 = time.process_time_ns()
    for _ in range(test_N):
        # Compute the (2^ea,2^ea)-isogeny
        Phi = EllipticProductIsogeny(ker_Phi, ea, strategy=strategy)
    print(
        f"Theta Model codomain took: {(time.process_time_ns() - t0) / (1_000_000 * test_N):.5f} ms"
    )

    t0 = time.process_time_ns()
    # Compute the time to push a point through the isogeny
    L1 = CouplePoint(EA(0), PB3)
    for _ in range(test_N):
        _ = Phi(L1)
    print(
        f"Theta Model image took: {(time.process_time_ns() - t0) / (1_000_000 * test_N):.5f} ms"
    )


def time_theta_sqrt(test_index, test_N=1):
    """
    Selects a kernel generating an isogeny between elliptic products from
    `isogeny_diamond.py` using `generate_splitting_kernel()` and times
    the average computation time of codomain and evaluation time for the
    isogeny chain in the theta model without additional torsion available.
    """
    (P1, Q1, P2, Q2), _ = generate_splitting_kernel(test_index)

    # Create kernel from CouplePoint data
    ker_Phi = (CouplePoint(4 * P1, 4 * P2), CouplePoint(4 * Q1, 4 * Q2))
    EA, EB = ker_Phi[0].curves()

    _, ea, eb, _, _ = DIAMONDS[test_index]
    PB3, _ = torsion_basis(EB, 3**eb)

    # Optimal strategy for the isogeny
    strategy = STRATEGY[ea - 2]

    # start profiler
    t0 = time.process_time_ns()
    for _ in range(test_N):
        # Compute the (2^ea,2^ea)-isogeny
        Phi = EllipticProductIsogenySqrt(ker_Phi, ea, strategy=strategy)
    print(
        f"Theta Model Sqrt codomain took: {(time.process_time_ns() - t0) / (1_000_000 * test_N):.5f} ms"
    )

    t0 = time.process_time_ns()
    # Compute the time to push a point through the isogeny
    L1 = CouplePoint(EA(0), PB3)
    for _ in range(test_N):
        _ = Phi(L1)
    print(
        f"Theta Model Sqrt image took: {(time.process_time_ns() - t0) / (1_000_000 * test_N):.5f} ms"
    )

def time_festa(test_N=100):
    """
    Takes data from the FESTA decryption routine and computes an isogeny between
    elliptic products using both the Theta and Mumford models, and computes the
    average computation time for the codmain and isogeny evaluation time."""

    # Finite field
    p = 0x176C11CF13E54B11406FCEC87BD4C1480F2BF6B3CF47C54370FEBD1C756E54F72C1501712922BAF5993402979D50DD13D09A841FED4773CFDB168F19A73E323F656921D7DCD797059B7B9AC3245C4D7BE6B343FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    F = GF(p**2, name="z2", modulus=[1, 0, 1])
    z2 = F.gen()

    # Chain length
    b = 632

    # Montgomery coefficients of the product E1 x E2
    A1 = (
        35894003914922960194491570668170134052925815028503174904348340271008987941056765900582123386036503409483517346111568603355977588827416898495013962955824223295844980918929024940478071024682704192582380112722991689881503452027891083904019566069064421282671925935011614895728854085690604796316303616857747672469672648234402242602400586496763226167634253173053471874825632413014160543913920750
        * z2
        + 23146170669210560502621292861089713471946652736458979067003229415535693206874706966935538994581754447944184064611672806470292750259663251209829011524312185911346092318222832612991175649417150715555908060248705993130270036556518615564774650228717302523894943958460593196903272851573247620550809134883751030874393916535285231402225449487865725181025901178245951836309471226776264689374198482
    )
    A2 = (
        33157640953836731568729867076929959670726332052161473833002218529005943652438605117445745658964481310158873521860644685031162573884926144401352383830514641134311115084066522185242808836860608411703598810503365396807631994717236543257811272497678811490006314722739948014914500940831668762316357121056673738616852633991320762551000069968411882299145268099642192114582770429942335243979487477
        * z2
        + 101997614125654320590390371380447161103016373745071029200573082749917015452092657758964741564090038665454586746959017255138531230278875005230861865341278061533003785291099337139896887799220097906387261529166855809739654566056674569446634688922834434195749365023939782814551722056758892066757569399321853491629860087517282943927135541762793575964076614796376361491889962198298512066240295064
    )
    # Elliptic curves in the form Ei : y^2 = x(x^2 + Aix + 1)
    E1 = EllipticCurve(F, [0, A1, 0, 1, 0])
    E2 = EllipticCurve(F, [0, A2, 0, 1, 0])

    # Points which generate an isogeny between elliptic products:
    # ker(Phi) = (P1, P2), (Q1, Q2) with Pi, Qi in Ei
    P1 = E1(
        33779709788973169224191220780030966979579809785320492444135379299185562317841147160011296942578004634705127528542475142703558480035569472565520162656647835955483546770450430235128274404484770618487878706273995263841527487549621804222670107544705589333887272598288291868099676263571499336190314295982868543898144571665886108500592228887659576632883901284152282591447104435353367069370609816
        * z2
        + 81684175990636284590332411665466061858566735003197010237165386831090729130925146937010533826699656419795420978842657443980139210299175677153978640361807604424415320816214351881778539832077326277155567133970547545227182495523058757913861332001470653847290390939960080238347804844039811783927698751642823183642237899610298671900451069540067537666365461698247800417647550647999869495971435379,
        71993418745746059209547489003349033302501317412583139435451250877703896123696599965656828554666009796467735881971669639060084306803231108463474758617314128379971834858681501845775930344303330076549459533245756073852088230543528241311533314735912832150385754804911711344208455913296915245990164905901045599520862129327499694971445723881029767922387078192629425986116105026392498570943451132
        * z2
        + 38760489756785639182787755884377364977030319409344438576499555759562122201727418515846220068918189414320002072287847774107159745448767024683674966874921396621843979911686905830315974257433383904422411726053799687017733346992447049348308713204787070476859634542339313762490710612213771321482592773978914461270563808831737556588742956428275961988539260412534746281719611465505536580702848860,
        1,
    )
    Q1 = E1(
        71349679068151173983737544530661483868442673835943419381098554738154613093749190510383613675715507213236722330634349749667121476547328943846038395180858425299538049454780196851498072085964113884727132309824443509646175567973527706878044535282245936090521852456166981425075023507551029392285045023151193725846384569196454554677587690204125883628842857249381570964083155718233463843956294059
        * z2
        + 106846856266064323120463401570949894666446363875754413212287254751084271089279144024274358464953407431553361262961439614920367207704389118198839341417369300335682816855898476177258739431383845519221998347401915975320600642323642103403807718865440330701692227998003595258194144462191843530299994569406005710788968961798926202233066357612663014238386723557799252125762428943310783023075077975,
        123586088799005214610093193595753262239666638790847322933344059369389774642476872851521527058505966869473263869385365437589397913046265845667977864314453183887873992977148912612119732611910638131573420971836686350916433232473305560594052355897309278918127100130093892417472650228453830318036895689421911938554416802729868422771295722227330220445171076934060905882898661196269706756879415319
        * z2
        + 10507140916115866437389730266182001053994248652785234259181280584714391587746712143097875420804449142861891251832953857210844276888951718780993600097052671714937910058014913619325749516366448615608599959873953917347455900817909636505538481307251790589644754095269575148711610137500288250029815137726638104673873865993682187535082905143317356573506051244678850409843282743985119644912273257,
        1,
    )
    P2 = E2(
        25013561878083992228194220941445311790806363908819025395085241260562042148587870970640510996611017390479630889340738510622072665692674503014258859306211095759752654528886054641740851025048926059441233413241170914485386166395409946032517774884102919125913608879776632151436336544200514746931395484653575908945635230083725126302809985311152643535395442167919997260823632003799878718176976720
        * z2
        + 83775823263801763058956181258264884582755536025908863521410991720480511856841366977660869014750689386110932303408278420448778965936736575492195242795428964957401079909352166812702849184400629569393268865993052400311147925799715658439788071020872982687139937128146105146879815880185082658373132213181487446393188070833815512888305841277489540391867544696604168946707256934457558798282555413,
        99824253066543114305562316003397858035144957108174019265782321203020208606855204672628457575516169535327187708197544661884552445587440896434077969213316747761481847462272894693681906889963126634963210058305152270748745447559016057099194465997247370519708917340465359543959407621827441914745628919812883035789821550877914776414992623757447199480379650973510691330955224783590189368709715251
        * z2
        + 78644052304901839233298904578832564342076494422883498859717679181325898411676130902959994463370686436619824204891733083767852186341600607390644716023395649153646645326187754110352286679838503184094338374867048991588832347197379583377335022849554210026991935131413017584167126277613453608714472625190724643795280601047552697606110664321871115469419132530117679457168472360753627209161737974,
        1,
    )
    Q2 = E2(
        113585246377665591945887763905879752875425348062392473701464618756825473340937984741547364993924990275392480726252607887492832851700567372402336366767958806853297541460712340734056992904141570575797551674839596701783161011641731513194117974893760867244265928991868830161057554113505465261148404250842022737775010403072998870937709036593551529441691382776086968812260131171502438614660462384
        * z2
        + 88464803224877467294939180751417989902006851603004076219567879912310113743662995530369558387691968541241703907648103277090968057590167330031086480083946243644205293664926712605685093434518785533005124902078683623325768253650524367422444196801229898474134669288984347782021956697380451320526967291777907865981321072821739367333353523749429995024227309801457600845303078109513529520163151176,
        111876406387148948193630898613296019618537750788309347069750540092350803571415174945605471389711449754121769773200009888955508948939687107318450721456536590010736500484541121673305794506245622880082965437018202661826016662812838610498603380570687869270931530213472906334097393022811316784234961350440353772310841611199328057479874204178756397716843852651172424126811604033984775799308288855
        * z2
        + 97307039093165336872797464570791375831784546965225917363736849569879456725064300315724332788662888564436692943928135050001759049687226845913183305770848611751425188364339983706059792692262083736036825819935813389337601499075937907371398456702144751047238995135010597526816554862117014843099164226665283256169153964469767038561814671898460274264920520050407134134772079272377899745431703623,
        1,
    )

    # Points L1_1 in E1 and L1_2 in E2 to push through the isogeny
    L1_1 = E1(
        75852364475916305027940748283088165885844105129068123624909260523486440602936600821992702429866923752891305381073929684630871248218079208891809821400327109921951004275009867529530487825180881962487099539980560338315936457366742512410466492280510053566407305524897232339240514921158880116994107209240733869102024018687443541434844962111866300932937989343678451473405035375050984253477138203
        * z2
        + 75683634740630257350636327250233438517263335330915694781687092496537368917395466398417746598220767447589516545081027719558717753385981050622192037129171362103181140373081395311203228376131859849147712727493983233768176403420045536656842350859231209772223690864485684910427745280296573852188538828733619820395047988130957748548113254776611415752130464538279080724888442256840345694129471064,
        65547734764658255763236604058660157023959404863677780021295213095558902871481751390864704389767530467075147763665962880597730690120872198624598314956746191976804941693757285389377433845184294714165048918451871913194903967436855709489214766391167857395738879639027491366098335235030835116931048796067293997326542964394200144329390956700744480369270867230258457207275708421426310796445378770
        * z2
        + 119631809608395499782676078115175356300895690588442037392934028116398442036937775159906182757656034212791669079153150154928577598521356616770094460488674584149635601467530040377333669590546595801377212289004403777521057250488781772762325148403514876391940594055658713543487385923659935060691078506275128157579407643855150116620787767706485535035641270437156239115174057794096815777589222593,
        1,
    )
    L1_2 = E2(
        91182349938370562123154916664394278362245126428711754887182239312713361158964721239517485193384428547918763880914691160861138482447858248699876472281754059244624814451358065926975806803743261859194233253896035669588395459839277737493866099235855532381323196326461899423422612180987497242381371861483039520346148521259132893767389508147789525593742544134373516881935367645028141919424099258
        * z2
        + 36255582990045705017272706125177967136462947478113400569429509717753674787915220650723925519755968612790683223346606015707330160910084300283727064975068326685087112620518481064797143475660567717323880314082150294062547524094035515867718461444121303582128973644008477823705666927622246463307158238981242035640183807992569406899021313569517886565762669422581837607928585145017384180104346966,
        25581700441851840972601547128885228966203797699824464699062502021644343918729943184213018363829820156766674256615137176168176895110648912153232162719002103212456178446220097027792574263743321481519521292371469548974283013362341328859304393631745764199627370305950873193431492589308456894553534574628599257327632187728526543987156916564929058978239752691456592702509399291901272999318076014
        * z2
        + 67449762652813345654638559973372090209648219693027334037342890783103831671691987673499876839443716857882067283701859812302452196545143030641838513195958322144872225430055793942086106772380281806718538687771220228095646788825019429900121505039719414056403628785431640695935475398571694188485199709256618941046117201462285111424153972049340561527636747785574888228821407057476339984891656733,
        1,
    )

    # Pack the kernel
    ker_Phi = (CouplePoint(P1, P2), CouplePoint(Q1, Q2))

    # Optimal strategy for the isogeny
    strategy = STRATEGY[b]

    # start profiler
    t0 = time.process_time_ns()
    for _ in range(test_N):
        # Compute the (2^ea,2^ea)-isogeny
        Phi = EllipticProductIsogeny(ker_Phi, b, strategy=strategy)
    print(
        f"Theta Model codomain took: {(time.process_time_ns() - t0) / (1_000_000 * test_N):.5f} ms"
    )

    t0 = time.process_time_ns()
    # Compute the time to push a point through the isogeny
    L1 = CouplePoint(L1_1, L1_2)
    for _ in range(test_N):
        _ = Phi(L1)
    print(
        f"Theta Model image took: {(time.process_time_ns() - t0) / (1_000_000 * test_N):.5f} ms"
    )

    # start profiler sqrt
    strategy_sqrt = optimised_strategy(b - 2)
    ker_Phi_sqrt = (CouplePoint(4 * P1, 4 * P2), CouplePoint(4 * Q1, 4 * Q2))
    t0 = time.process_time_ns()
    for _ in range(test_N):
        # Compute the (2^ea,2^ea)-isogeny
        PhiSqrt = EllipticProductIsogenySqrt(ker_Phi_sqrt, b, strategy=strategy_sqrt)
    print(
        f"Theta Model Sqrt codomain took: {(time.process_time_ns() - t0) / (1_000_000 * test_N):.5f} ms"
    )

    t0 = time.process_time_ns()
    # Compute the time to push a point through the isogeny
    L1 = CouplePoint(L1_1, L1_2)
    for _ in range(test_N):
        _ = PhiSqrt(L1)
    print(
        f"Theta Model Sqrt image took: {(time.process_time_ns() - t0) / (1_000_000 * test_N):.5f} ms"
    )


STRATEGY = {
    126: {'flag': (True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True), 'doubles': [125, 48, 30, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1]},
    124: {'flag': (True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True), 'doubles': [123, 47, 29, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1]},
    208: {'flag': (True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True), 'doubles': [207, 79, 52, 29, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 23, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 9, 5, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1, 1, 1, 1, 30, 20, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1]},
    206: {'flag': (True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True), 'doubles': [205, 78, 51, 29, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 22, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 8, 5, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 30, 19, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1]},
    632: {'flag': (True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True), 'doubles': [631, 265, 152, 86, 49, 32, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 14, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 5, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 66, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 113, 66, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 47, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1]},
    630: {'flag': (True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True), 'doubles': [629, 265, 152, 86, 49, 30, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 66, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 113, 66, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 47, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1]},
}

if __name__ == "__main__":
    print(f"Testing 254 bit prime, with a chain of length 126")
    time_theta(2, test_N=100)
    time_theta_sqrt(2, test_N=100)
    print()

    print(f"Testing 381 bit prime, with a chain of length 208")
    time_theta(4, test_N=100)
    time_theta_sqrt(4, test_N=100)
    print()

    print(f"Testing 1293 bit prime, with a chain of length 632")
    time_festa(test_N=100)
