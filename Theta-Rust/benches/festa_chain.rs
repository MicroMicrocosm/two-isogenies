// Benchmark the performance of the isogeny chain of length 632 between elliptic
// products with a base field Fp^2 with characteristic of 1293 bits. This is the
// current (2^n,2^n)-isogeny used in FESTA-128 during the decryption.

#![allow(non_snake_case)]

use criterion::{criterion_group, criterion_main, Criterion};
use std::time::Duration;
use theta_rs::thetaFESTA::{product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, Point};

fn criterion_product_chain(c: &mut Criterion) {
    // Curve coefficients
    let A1_str = "d268d78b4a92566637452bf184ade847911397b523d84750d9beddf94d402fa5ca047b8a729844f8362094d212cd6fb9121a44241fbc8a0cc2bac1562744f015e3d41a6fa75014fa6a00ad4ac9ea34b202160d0876574324a65b7b2253142a6f408e89094e8ed1ba44fb56eae37cd6de728f177c83391c6763bbf493fc33874cee9259989d59d403886d48dad157f40238acb776ecf4e38fb5653ddada4afd6d0402099fc10c379b5190c8c195531815b4729da27d43cbaf41b5d53185cabaff515669eb38ef7a39eb1ecbae997f063f47fa9952aef62eed5ae8b647d578dcc1da25aaf197d1c2cfcb9043e9751c71fcf42ff71d5e975279fe6d0d1ee051b167d5e23b86f56b4cf99b5c52a8b22769b9b7c3d7d224afb3cb3e77d89837765a983c5cddc1e3174f55fbe955a3a4b13d87d9746edce0c31d1bac44eeca2cf0973c240c5b0c";
    let A2_str = "3d4cfe9cbfe6987fb8ea7a04f62c9ed28ef9f8ced54ed453efcdaa3621d00b77b0e7a7dc12b1d582032f0b98e174aa8dcbc514b88f39988a54e6ae6be1cf4f8746a86a8837fa6f10fb999732616484422cdd7ea7e501b858f6ac3f4cfdcfcd467301ec08817a9a2955e274391d16f57d0667eadbbbe77c7b48509ca21b5ccb254d879f13b10a2dd3d089aea44ecc8bb5efac9764f86630460a5a9ea718b4394dbb091ffa5bc28624f6004c4068a7b9dce26d84d02a5a0b50b1bc168d7b9ab63b920a1ab0f7b35c09e6ed717331281699648e4f723916033dc3784491c12470bec626c3a689b95b746c342cae065dd8a8592949148d0a4a38ce96238cb6c6a5bee6e272e67462a0afa95c4f06ad7151ae1742483d1ad3203718650d40907c717eed90dd7b4e3479b623b973cc1c981f57ea2c24e76a13c982a184a4252cc3e80d0d5a5e15";
    let (A1, _) = Fq::decode(&hex::decode(A1_str).unwrap());
    let (A2, _) = Fq::decode(&hex::decode(A2_str).unwrap());

    // Kernel Data
    let P1_X_str = "ab710145f609f21b96a9829e5e9ca411431f13fe7155e1a6424013967a326a691c701d3fe90cf1536cf504fa70e0c49fcc401b7e104bec8f5e648b64bf54694c83377ac40e795be89f394e2c581a79db88f9d58d5c3114d4c07488ed95467747ac88add4ddc2f8aa95eb47be54c115b294b34fa8fb8858224b5289cccb1b1fbef87cbb7b9f1acb012c9e4341b40e7db00bdb4ffcfe88a27175aa488e887c1cf1880cdaa0a1853d95790a0902169eee111a7cf76765d92f803f615a3fe73fc01d9d245bb0325dcb5f13017c344206b583d3324fd466509308f329c3e40a6ed152e01d763eeebc2525a3e1056f82831111b54684ebcf9776d4f0be11a92de06b5c47fca1188579177bca071a152e16133e17a055da98e6bc706557f35942c0e3f665a2db4c489ac91cbe484b3092fffba9b545799fcd59626fbfe5fb77c47aef6aeb5b7616";
    let P1_Y_str = "c05a766bc0b28b4977118f646c78d43a26c1f8b7ae75847f63f6ed94e4247d581e65e85a4d86eebbe8f7928670e8962e3c51cd843b7aec237c5e1d1a9a06e7a74c7a3672e0fd8f12704bb95ec930599e5dbb2790c9480f345feb74acb27050427107aba7227730e0a6f1951a3ff22e459ef26786d93778d6f464c47ff85773985726614a91f5da610a3cf5aad29f59bd70a7c3e41861eaa0e3636df16ec987f27816753310e1ad9014a3887b2e86055cc0b087fe71148fbef4ed961df36f7a74920904709290ac3215f35056257fb284c08b02516c6ba3f7c5b66001055f0f606fad2c24b4ab975012b10ef589e1600e90f02decefb45895aa28f70d931ccce1cb17a7a5383189e10927ec2ab1dc3fb159496c7f947fcd8f6159a04ae8fd1d401929978ebb0e7d2a7ea554cfb9230351d3b39f9da0975477e17ce9291f807fc39be24813";
    let P2_X_str = "495a0140eb8101aa90f13bdab5c65e53a655ef2b8d2280b14396f498dad675305b822928d54d6f81d494d5e639e28d6c29585e788454609dc69653e5c1294b6ffa69c17e5737a8bf506c0da45a792140d167e8d67fe0b39288f0cd67cb33989f147a8d88189967e384eb54d40b105424eedd368b89b8d26beac0ebbdb2fc2bb12aa105116ecd2d076fdea288f15af4d89f543a85077db767e2ca19c7c94f3331330d5d888cd541b6297558e74cb8da409ffd041094b0e8f6725fdcda3052dae881bf1f53bfab770946cc8167ae8dbd0c5352c4153a64f2e9999bd31ae1a8ce27ad9816d5dea6bb7101a856c3804e3d26eee3b5656c58ba43a00460a9fbd4f808dd21e998abe955513a11ee68e2bcd80769ed5a1876e2f7ddac97b50d790b494abd1aa733ff76956ed56fe2457c9175a8504a52bbc9ce1ac21b9a2c2b6fc517af15b59307";
    let P2_Y_str = "eeb2858418f1768adf8e61b803381587a72dd971a982a2cf40cc263905feddf6fd4131a134c89bcc46f8a758fb1ec03dfe4d973fbac5f2533774adbe8ef0be3b7efc4650923e36f79d906d443ec7fef2e97e1dd92189ebc0e2a29556dfa821508b353b48e8b231785ba7c95e38b9ad847b08e2949be98eacfdfaa3f7837b5c436dc3034d78b6d975bee83943981c25cabe5f33a6a39faea86c4501989cd51516f6033be07b5580b03020820669a0d84ab8456837f90c665526b41ec7fe76d3a4e0c8133890845b5b92532cfdd76c1171b2b824a9bbd7b5e1895798d313f2145abd944b91d7a0d4407b608584e72ac91e4d24472f670430bc7c7bf7c6e9829a4063a254d6dd245c52fe833494b430e4cc01f34a50d823a908427908f7ddd7704bbd0ad34bb7e3ae1d10061d214d6020ad2390730393ce8fd9da5afa8252923b88463ce20e";
    let Q1_X_str = "f1327865b3ecd10fd9af8a10a7301c031637d39d66da8caa9952cf5b9e1ab8b518303a2eb8dc913bef7fb352242e3d80b27056eedded2bc1fcf5e70bb5afd610ca066c1e34c64d8801d98af75c459f805a03719a7ac29ce07b5e5625a8eabab6c705b87d217fadee068dec01c9adf4071f5eace53cc7a50914af62dff0d5343c5efb28f67fa3f4018a7192eeb0ba11bcc034f7425bce1adb07369c92268a1f3d170547805c63ef9f44f099cf3f1733304fd012a3feec1ba33c7f9626e24c83ca089ff82005e39724dae7202ea80c5e3d1271b2aeceb6046dd26d01228385f8fb4493c05d0c0a6f9c3eb9fbe4ed0cd7dd2e9bf0896e6604b4a2f9e3b4ce8b7eb8330002f497bba8985854b8652b29a55b41af08ccc4e3db74d71d65ed7ecdef8ed821ae3ac65bfa60f23b089b674d21413a6bfc04cf0345ad6ee50d16a65b84e83564df0c";
    let Q1_Y_str = "6c3b1385160d8c41e4fcd4b0f12908a0ad58772a35b34fc0a15159a5a164ebca6f53873d8b44ea79435e1d3c7afa318286e2c952070f81cb46d8e55a9ce221fc0637e5f3e5dba2422e4858e3ce6720f5070e7ff93837c2d1df96327865497dfde7160e7e6fce952ddf6a4c355d0c65127840d73999e7b81a447d07e694d873639957fb67c8d3f5e10622a03e151899ec68f24c1879fa8f9e44f4749014e27990e60272fa04158db4ad5b7274eb88a456ed34d5e67cf1eaa8867f2158f85355c5ec68bb04b086622d28ef2ff62dfde9988d79e5444edadb59fbc999c9fdaa33488c94b178d9c218109682ee4ba51b8709fc81bb66f0683e219b6c5b9f2a7f05dda514329645d7858fa65093a86fe8864bf1dcabe9f5d3c952c27ab45a4491d9291f37b2876bb62e9d96ab3fa29bb5133c18ffab247e024cd30b23f2111eb7d4a3a9759b0a";
    let Q2_X_str = "34678391557c9e2776dff7177e4c6fd8e62442ecc0aceea3a2b7969b9720c493fb4d96affe30264e5ab4c5e73f63d35da8ccde60fe738470c54b849fa9a99d78a415e651c583f18990ce9483cc905c2338dd704b738c9158febea12db9bf2c546a101a0fbf6e4923f9fd89090ef4843a2005dab3c0a14f080b7589cec37f7ee9477da4cdadf3311ae9a599734ff2de0e9fb364627eb0f5ce10ec05cadd9f8e4fa904b646f943c1cd52a1ec46145364f1edc41700525db027b0014b97caab13f747fa72552e51e41e185064178d643ce49547eff8c06bb232ab372f6d4ecfeb8205ca7c42f70f0e8085d15a91d3801bdec8aee5ea1c516410df8fd3f0f09be42876614fc25ee675fc743d1d7f58b6a31af37370a5002ec639ed014f4a2c6fbf1b646b6a24b6cec49e3c9f7e1e1095e32f880640f47ed835c50aa7ea313d9068f1d7e6980a";
    let Q2_Y_str = "7c33493b0c0ab5bf88131c9bc9291b428d6812683d1f82b0f6c41dbc95e999b336a0a9efd0050255f0f290d76291ea7c72bb03d3e3130a2101ea5b4dea3a62397f56242e5ace05bf2754db2bd51ad90f37d92a4b7fdf9e705b4788b2a9e6088faf446475bfeb4cc8fd19b6655c5953611e673d3cd354132cf484e9e395eb557991ddd4d34a63036ac95801769d2ded726c8c77ae8d9bdc29dc42c34a75d03e0489047870e6f93ebcb5ec4cdf849b2c7f28404f3380c9b31b02e43e3c8a1c5e892aaa68d1c260f0d903055021027db1ac829ea2603f775e9c3b6d5594accb05c0656d9dcbd9c030510d1c3de017cce56e16d16a4e8a0aa5f31e1d9fc565cedd7757bac076f4d06fd9b24887af3d5829a0df783f0485890f9c1c1446cb99dedd1ddb81beb6f48f8c6a1928aa053d9c3c634ea104388d28a592e4c2aab38fd01999aa797c0b";
    let (P1_X, _) = Fq::decode(&hex::decode(P1_X_str).unwrap());
    let (P1_Y, _) = Fq::decode(&hex::decode(P1_Y_str).unwrap());
    let (P2_X, _) = Fq::decode(&hex::decode(P2_X_str).unwrap());
    let (P2_Y, _) = Fq::decode(&hex::decode(P2_Y_str).unwrap());
    let (Q1_X, _) = Fq::decode(&hex::decode(Q1_X_str).unwrap());
    let (Q1_Y, _) = Fq::decode(&hex::decode(Q1_Y_str).unwrap());
    let (Q2_X, _) = Fq::decode(&hex::decode(Q2_X_str).unwrap());
    let (Q2_Y, _) = Fq::decode(&hex::decode(Q2_Y_str).unwrap());

    // Points to push through isogeny
    let L11_X_str = "dea2bb5f803a0281b2c1dc948ce739aa77ccb66d167b62b069c4e6c5e8f2d6a21bb7b385d134f04608caf3f129aa29ed6ca8be9369ed27c43535bd282878c1e8c8907581110a22df7b2ab1358b8aa50fc665b7ca399aea781ef0ac250277583cb97dc7116f98e6a45810c468a9d5da4825f3af65c4a6b96debe9fa0b5f7805312eb5403a7fa97da650e475d1daa6a93db1e5bfef4ef4e977fae52f416a2af1cf2c04810a0c60deebb39bf19d01e74e0ed02cb85518d7285d804b56340e7e68bd3e59312c2b4cad5ece56415b5ec824c4be694df9255ff04cb3905af56adf4d6473efda1922ee615d1d3c0741eadf41d1c9587d49ceccaed909b8b0174e61ce6f0df628a03da0def573fd3b7d4850bfa751b73a37839ebb484f9a9f112d95d1ee2522dc3703c5102d957fc3dfbcc24d5618a6cf5b106ef1f23590fc603699aa5887f01302";
    let L11_Y_str = "abe10e6a88932185ab71f34c0a808b1a6e34844f6bea614c04515b2755f84947d13eb29658659a4d4c66a0e2031964448f69b14a3cfe8a08533981747a82d1e648bfd9faba0b5349fe1a057e356b4fc9fa6d71078173665cfb07e19b03e2da3c45aab7e70ca1c80a0a9e0bcb2b551103c58da65277897c2470b70a05595ade7f7950c04f57497c8421d9081591529f54942fd749fb86c28ebedfd85e5d3e8f92131275c689849076add344565b166beee4d889791d9072f313305ce5c8c1f39c6c339e986d852672a1eff6a7f08ed7c99abf884db79806adf520c2fc8a34d6ddffa600454f181b563375de64f62e16ac48aa413301d5d4443ad325552d3381c4b18e317ed659978967a71d3a68c9238b6828d4f9483f5c298a872203036e5380112b1da86f0f7b44c205d20c1290eec4c644f0b5aefa1dbaac61142bc8d4232ac10cc30c";
    let (L11_X, _) = Fq::decode(&hex::decode(L11_X_str).unwrap());
    let (L11_Y, _) = Fq::decode(&hex::decode(L11_Y_str).unwrap());

    let L12_X_str = "67bd25d443ab35439cde45a4cae9deede7c7ab808d3cf080ae2d2fa8ddfbbde33fb472b24d0838a80ee4af7f4d3fe5ae85a1becad033785fb7024b4d203b3372c9a71816f0f16108fe4f9f989dec5d8427af9f49a7f3c43c3f7e8eef91bd8126cdaf3e953ae90fcb679e9d84beca4b803af9f8f6ab0af82ab00693a258121450afa5fb59227855a0e6b7301362382b8bd040d72d5e48a5627bd73d1b8c4c3e94d00fe3807cf8708a8e5752c0c4b42e9b68cf6640bda42a05efb996560aac1c407a6c0b36bbbcb9c6355cb897bf7191aafdd8ab7582becd10a5f57e40851debd9076ad86d4c425152e8a5609ee2a2c50737d05e09a2089cec0066d6e8a8f8e613aa06cd268e09725ca0ea0ad865753779e2350f7324a41145f57eec8ca3527fcd28024cab1cf831ae25398f88ad3052cd6a8955573af4e79e77eea6f3353282234072140b";
    let L12_Y_str = "fc5fb5ec8552da8cdde3fb67034abcd9832e2ae0ee73a46de6d17d1f3565eeaffb053011a34e951572ac8ba83840f713c104f171fd76979ccdb3fdba531df2fe172d37bd01a7246a3132887c6ec5cc76a6e8a031df962a86bbbdfc0ec813b501a1246a34d784a4c4bca9cb10785a7eea9a068df98907769c351292287cfa6f441de196f2406e966ab412140789260f8c918390c8149df0513bf3fda0eb1d275dcc00d30f507ce8ad17a91cdb083227b971162b04284eee36506ccf6c71b1ee07d3f2ac06e74192cd560e3d6093a4375b18b7b6281c9e1cc0583495e5169e85858ac2be3e600c6a8ba53f6c7bd9ede931fca251fa0fe2a218ae3d47565d45fdcdbaea8a262a22f9b472ee2c6fb3df9963fa2a6d8899485659ea42110f0b5dc8c92855de8dee7e64ebae6c45af01f83e36a7a3f8775808867cb165e04918525889d1661706";
    let (L12_X, _) = Fq::decode(&hex::decode(L12_X_str).unwrap());
    let (L12_Y, _) = Fq::decode(&hex::decode(L12_Y_str).unwrap());

    let L21_X_str = "5b3bcf51421852af967115e3824bbeffe9b9aab16762e530e35c32a77b0b37c62a7e0ce96f7af6b831cf7193f807ed48749ebee12c5fcf8a2fd9bf187b5a58856888f99c81fb89644138d0b136589aa28da0feb085cf3d2a43c0d54bd8c87aec9b818542ff96c0a95643402026f29eba75c9cdf0deabd841bc1f1e88f5dd04bdac64174c5b7dc7af8428f3c74c828d5dd37ff1afda9e00d4a4b8dc3d579bb1425809d8a751619517a25f9ac408ddc3391c73b170127d1e2492d8877ad069052b6cba9159ff2590bc8e628496002bca3665174d414ed67d842f93931418560264ebb0359f9ae42ca64ddd75834cdeae6e7e7172f981a419406fe24401c4d70198cfadb6f301125b7c3f37a65babab3f347c935d0cbf8ef8e7dda642767960e73380fef72976e386007bba68fefc04ab876a50db3e0c7f9b12965b61e4fbde8f887c47f70d";
    let L21_Y_str = "5ff5ef87185ee6951c538ef7f5e16b17c3c9a28b8e7828743e45472482d01062d5d616cfb34d86fb0a66cf65aa7820353b7d8df025e91b6685009ef4b2231983c3b428e79e11afe00aa7f69a7b4f4f1c09f51953b44ba84b2b0d5e80c47f66a5e883353acca8b79cd65b5c77185c028a5af8b08cdc33b01ae35be1866b932cce8ef2803800a50f07dd6d835503102ed21871f13814051ae73d8d1e29ef85a94032082996a0c8e58c5460b9dfef65311511e8030bb434ce7e4bd98ee7981143f2ac175e141776e583ce6ad2e2e06bc6ad2e20effe7969b9cd57ac1adb7cd89decce7347fcf1522dfed1ff2663d3942fe9fd86dd798d9ff6edb0db17472259e5b63086e4e8e04f112c474f8a4b64dd5a89cfb176ea9742f954a507b21dc3163bc1609dd5077daec60c14851496d3e56612210a8bff825890752462c02e8b720d909b685308";
    let (L21_X, _) = Fq::decode(&hex::decode(L21_X_str).unwrap());
    let (L21_Y, _) = Fq::decode(&hex::decode(L21_Y_str).unwrap());

    let L22_X_str = "4eb3160cce9bb8a02fbd87d6f070cfdab2f0465f2ed692457f56a9699f3d2a51385d1b0ea6c112676da50db2cd916295aeae5a0c8fab327c9019ecdbbae67d0677c38688b7988a99bd07bc4f67f47de987f535bdf4f7af89b6f17aebd5841eebc675a3bd0622b4445a86f469ab6478760e8d54f4810e60f82f25df97957768708b2373e28c5bdc051f5eb829ff63e25f024f546d7c0d2d30a90d235961ff861f14163fe43d8b8c9b8fc5f8bc8ee7970f3aa23ab78eb4a9a836ff087806843a4935860aab038b0809f7df120f5eeaa5fcf6c936158448cac458ac9c7ca3e89d830efde790240a79a515c4cb76e6e93ed2b121bbb75397153154ed16345358a9b2e6e28c70cef892dc01fb13ee017d04ee3848cf67fc8107b80b02b139ac9cb987c47cc3b8746a49d9e21e68643861d94cfbdce8aaed2ec3fd81be6b4aeefd27a3d85d7c12";
    let L22_Y_str = "6eb74be4d694cfd39f3a216ea77514079c66b8dff598af6588683946493a691d812685e352c218f95f465e2b21da95012c409cbf5c9925ae5e2c3ba6249164641d0314f03de05f3f9cd6851ae2681f04a7d4d04a76999ff6ce0dcb73f2d1e76fe9c86a403c62f4de281c3bccbfa72410a9bf26a7d2959577be8d9b23de9c634fd096b468e1da9ee98a4a449ac99bd590c1f9b60615dc748b6d7966e61e373edb7b07c97cffcd09ebc8daa380977a10555f800229e2bd158acc311e4e0bda66dae5f2b6d0745d62ed3d9a7d5415ba6aecfcc7db6eeff54476179c59e5cdb6f62c32f634b2de573181629b25d3febea4943982f72b535389183118caf16685fdb5ccfa73e609f9a7a98b7b1209231d0d9004cdff6b425cf578a4f586473e9a4fbf32b6fbae0cd1ecc1ea48d6ec950b2c42fe9a841a9cd559ab0154dc60b73ea290d1d37613";
    let (L22_X, _) = Fq::decode(&hex::decode(L22_X_str).unwrap());
    let (L22_Y, _) = Fq::decode(&hex::decode(L22_Y_str).unwrap());

    // Curves which define elliptic product
    let E1 = Curve::new(&A1);
    let E2 = Curve::new(&A2);
    let E1E2 = EllipticProduct::new(&E1, &E2);

    // Kernel Points on E1 x E2
    let P1 = Point::new_xy(&P1_X, &P1_Y);
    let P2 = Point::new_xy(&P2_X, &P2_Y);
    let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
    let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
    let P1P2 = CouplePoint::new(&P1, &P2);
    let Q1Q2 = CouplePoint::new(&Q1, &Q2);

    // Points to push through isogeny
    let L11 = Point::new_xy(&L11_X, &L11_Y);
    let L12 = Point::new_xy(&L12_X, &L12_Y);
    let L21 = Point::new_xy(&L21_X, &L21_Y);
    let L22 = Point::new_xy(&L22_X, &L22_Y);
    let L1 = CouplePoint::new(&L11, &L12);
    let L2 = CouplePoint::new(&L21, &L22);

    // let image_points = [L1, L2];
    let image_points = [];

    // Length of isogeny chain
    let n = 632;

    // Precomputed from strategy.py
    let strategy: [usize; 631] = [
        631, 265, 152, 86, 49, 32, 18, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3,
        1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 14, 7, 4, 3, 1,
        1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 5, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 19, 12, 7, 4,
        3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1,
        1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3,
        1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1,
        1, 1, 1, 1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1,
        1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 66, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1,
        1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1,
        1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1,
        1, 3, 1, 1, 1, 1, 1, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3,
        1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1,
        1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 113, 66, 37, 19, 12, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1,
        3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3,
        1, 1, 1, 1, 1, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1,
        1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1,
        1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1,
        3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 47, 29, 18, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1,
        1, 1, 3, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 11, 7, 4, 3, 1,
        1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 18, 11, 7, 4, 3, 1, 1,
        1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 4, 3, 1, 1, 1, 1, 1, 1, 1, 1, 7, 4, 3, 1, 1, 1, 1, 1,
        1, 1, 1, 3, 1, 1, 1, 1, 1,
    ];
    //let flag: [bool; 632] = [
    //    true, false, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    false, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, false, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, false, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, false, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, false, true, true, true, true, true, true,
    //    true, true, true, true, true, false, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, false, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, false, true, true,
    //    true, true, true, true, true, true, true, true, true, false, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, false, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, false, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, false, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, false, true, true, true, true, true, true,
    //    true, true, true, true, true, false, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, false, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, false, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, false, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, false, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, false, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, false, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, false, true, true, true, true,
    //    true, true, true, true, true, true, true, false, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, false, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    false, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, false, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, false, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, false, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, false, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, false, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, false, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    false, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, false,
    //    true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    //    true, true,
    //];
    let flag: [bool; 632] = [false; 632];

    c.bench_function(
        "Product isogeny of 631 steps using optimised strategy with 1293 bit FESTA prime",
        |b| b.iter(|| product_isogeny(&E1E2, &P1P2, &Q1Q2, &image_points, n, &strategy, &flag)),
    );
}

criterion_group! {
    name = benches;
    config = Criterion::default().measurement_time(Duration::from_secs(15));
    targets = criterion_product_chain
}
criterion_main!(benches);
