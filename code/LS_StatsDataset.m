
function [anals, toymodelfits, preCalc] = LS_StatsDataset( config, fdata, focals, toymodelfits, calculate ) %, fcollated

% temp = [];

if nargin < 5 | isempty(calculate)
  calculate = 1;
end




%% calculate simple genomewide average statistics from data, if not supplied (these statistics do not depend on the inference)
iC = config.chromosomes;

anals.gwa = genomewideStatisticsLight(fdata(iC)); %genomewideStatistics(fdata(iC), 0)
CSChrLen = cumsum(anals.gwa.chr_len);

for c=iC
  preCalc{1}{c}        = calcPreLHStatsSw( {[]}, [], {[]}, fdata{c}, [], config.predef_idx_test{c} );
  preCalc{2}{c}        = calcPreLHStatsSw( {[]}, [], {[]}, fdata{c}, [], config.inf.predef_idx_train{c} );
  %   anals.gwa.EgMutDiv(c) = preCalc{1}{c}.EgMutDiv;
end

for set=1:2
  sSamples = 0; sHet = 0; sHom = 0; sSDiv = 0; sDiv = 0; sSamples_hc = 0; sHet_hc = 0; sHom_hc = 0; sSDiv_hc = 0;
  for c=iC
    
    grat  = spatAverageRecRat( fdata{c}.genmap, preCalc{set}{c}.pos, config.rec_spat_window );
    
    sSamples = sSamples       + length(preCalc{set}{c}.idxc);
    sHet     = sHet           + sum(preCalc{set}{c}.samplesHet); %preCalc{set}{c}.samplesHetSum;
    sHom     = sHom           + sum(preCalc{set}{c}.samplesHom); %preCalc{set}{c}.samplesHomSum;
    sSDiv    = sSDiv          + sum(preCalc{set}{c}.gMutDiv);
    %   sDiv     = sDiv           + length(preCalc{set}{c}.idxc1d);
    
    sSamples_hc = sSamples_hc + sum( grat >= config.collate.rec_pol_th_L & grat <= config.collate.rec_pol_th_H );
    sHet_hc = sHet_hc         + sum(preCalc{set}{c}.samplesHet( grat >= config.collate.rec_pol_th_L & grat <= config.collate.rec_pol_th_H ));
    sHom_hc = sHom_hc         + sum(preCalc{set}{c}.samplesHom( grat >= config.collate.rec_pol_th_L & grat <= config.collate.rec_pol_th_H ));
    sSDiv_hc= sSDiv_hc        + sum(preCalc{set}{c}.gMutDiv(    grat >= config.collate.rec_pol_th_L & grat <= config.collate.rec_pol_th_H ));
    
  end
  if set==1
    anals.gwa.gwP   = sHet/(sHet+sHom);
    %     anals.gwa.gwD   = sDiv/sSamples;
    anals.gwa.gwSD  = sSDiv/sSamples;
    %     anals.gwa.gwP2D = anals.gwa.gwP / anals.gwa.gwD;
    anals.gwa.gwP2SD= anals.gwa.gwP / anals.gwa.gwSD;
    
    anals.gwa.gwP_hc   = sHet_hc/(sHet_hc+sHom_hc);
    anals.gwa.gwSD_hc  = sSDiv_hc/sSamples_hc;
    anals.gwa.gwP2SD_hc= anals.gwa.gwP_hc / anals.gwa.gwSD_hc;
    
  else
    anals.gwa.train.gwP   = sHet/(sHet+sHom);
    %     anals.gwa.train.gwD   = sDiv/sSamples;
    anals.gwa.train.gwSD  = sSDiv/sSamples;
    %     anals.gwa.train.gwP2D = anals.gwa.train.gwP / anals.gwa.train.gwD;
    anals.gwa.train.gwP2SD= anals.gwa.train.gwP / anals.gwa.train.gwSD;
  end
end


% return;


%% collect diversity, divergence and properties values
poly_gx=[]; poly_px=[]; poly=[]; DN_gx=[]; DN2_gx=[]; DN_px=[]; div=[]; sdiv=[]; recrate=[]; recrate1=[]; Dn=[]; Dn2=[]; recrate_Dn=[]; ocodons=[]; codons_gx=[]; allsites_gx=[]; codons2_gx=[]; allsites2_gx=[];
poly_gx_shift=0; poly_px_shift=0; DN_gx_shift=0; DN_px_shift=0;

for c=iC
  if ~isempty(preCalc{1}{c}.idxc)
    
    polovtzi = zeros([length(preCalc{1}{c}.idxc) 1]);
    polovtzi(preCalc{1}{c}.idxc1p) = double(preCalc{1}{c}.samplesHet(preCalc{1}{c}.idxc1p))./double(preCalc{1}{c}.samplesHet(preCalc{1}{c}.idxc1p)+preCalc{1}{c}.samplesHom(preCalc{1}{c}.idxc1p));
    %     divmask = zeros([length(preCalc{1}{c}.idxc) 1]);
    %     divmask(preCalc{1}{c}.idxc1d) = 1;
    if isfield(fdata{c}, 'Div')
    DNmask = zeros([length(fdata{c}.Div.isdat) 1]);
    DNmask(fdata{c}.Div.idxN) = 1; % [fdata{c}.Div.idxN; fdata{c}.Div.idxNfake]
    DNmask = DNmask(fdata{c}.Div.isdat);
    end
    
    cur_gpos = applyGenmap2pos( fdata{c}.genmap, [1:fdata{c}.chr_len] );
    crr1 = spatAverageRecRat( fdata{c}.genmap, fdata{c}.Poly.pos(preCalc{1}{c}.idxc), config.rec_spat_window );
    
    gpos_min = fdata{c}.genmaplims(1);
    ppos_min = 0;
    
    poly_gx = [poly_gx; poly_gx_shift - gpos_min + fdata{c}.Poly.gpos(preCalc{1}{c}.idxc)];
    poly_px = [poly_px; poly_px_shift - ppos_min + double(fdata{c}.Poly.pos( preCalc{1}{c}.idxc))];
    allsites_gx=[allsites_gx; poly_gx_shift- gpos_min + cur_gpos' ];
    
    poly     = [poly;  polovtzi]; %polmask
    %     div      = [div;   divmask];
    sdiv     = [sdiv;  preCalc{1}{c}.gMutDiv];
    %     recrate  = [recrate;  fdata{c}.Poly.grat(preCalc{1}{c}.idxc)];
    recrate1 = [recrate1; crr1];
    
    
    
    if isfield(fdata{c}, 'Div')
    DN_gx   = [DN_gx;   poly_gx_shift - gpos_min + fdata{c}.Div.gpos(fdata{c}.Div.isdat)];%DN_gx_shift
    codons_gx=[codons_gx;  poly_gx_shift- gpos_min + fdata{c}.Div.gpos];
    %     recrate_Dn=[recrate_Dn;  fdata{c}.Div.grat(fdata{c}.Div.isdat)];
    Dn       = [Dn;  DNmask];
    %     ocodons  = [ocodons; fdata{c}.Div.isdat];
    end
    
    % fusing genetic maps of consecutive chromosomes, with small break intervals
    poly_gx_shift  = poly_gx_shift + diff(fdata{c}.genmaplims) + config.gChrSpacer; %0.1;
    poly_px_shift  = poly_px_shift + fdata{c}.chr_len + config.pChrSpacer; %10^7
    
  end
end


% if config.smoothgenmap
recrate = recrate1;
%   [sY, nX, comX, eX] = poolYbyX( poly_gx, recrate, (max(poly_gx)-min(poly_gx))/config.smoothgenmap_win);
%   idx = find(~isnan(comX)&~isinf(comX));
%   Srecrate                     = interp1(comX(idx), sY(idx)./nX(idx), poly_gx);
% %   recrate(~isnan(Srecrate)) = Srecrate(~isnan(Srecrate));
% end


if isfield(fdata{c}, 'Div')
% config.sDNwin      = 3*10^-4; %20Kb
% config.sDNwin2     = 3*10^-5; %2Kb
% config.sCODwin     = 3*10^-4; %20Kb

% TARGET:   SCoding, SCodingDensity
% coding sites in cM  (sCODwin)
[~,~,~,~,~,SCoding        ] = smooth01NonUniform( codons_gx,   ones(size(codons_gx)),   poly_gx, config.sCODwin);
SCoding                     = 3*SCoding         / config.sCODwin /100;
% sites in cM  (sCODwin)
[~,~,~,~,~,SSites         ] = smooth01NonUniform( allsites_gx, ones(size(allsites_gx)), poly_gx, config.sCODwin);
SSites                      = SSites            / config.sCODwin /100;
% density of coding sites,  sCODwin
SCodingDensity              = SCoding          ./ SSites;

% TARGETS:   SDn, SDnA, SDnB,  SdN
% probability for NS substitution at a CODING SITE  (sDNwin, sDNwin*2, sDNwin/2)
SdNcoding                   = smooth01NonUniform( DN_gx,       double(Dn),              poly_gx, config.sDNwin  ) / 3;
SdNcodingA                  = smooth01NonUniform( DN_gx,       double(Dn),              poly_gx, config.sDNwin*2) / 3;
SdNcodingB                  = smooth01NonUniform( DN_gx,       double(Dn),              poly_gx, config.sDNwin/2) / 3;
% coding sites in cM  (sDNwin, sDNwin*2, sDNwin/2)
[~,~,~,~,~,SCoding_forDN ]  = smooth01NonUniform( codons_gx,   ones(size(codons_gx)),   poly_gx, config.sDNwin   );
[~,~,~,~,~,SCoding_forDNA]  = smooth01NonUniform( codons_gx,   ones(size(codons_gx)),   poly_gx, config.sDNwin*2 );
[~,~,~,~,~,SCoding_forDNB]  = smooth01NonUniform( codons_gx,   ones(size(codons_gx)),   poly_gx, config.sDNwin/2 );
SCoding_forDN               = 3*SCoding_forDN   /  config.sDNwin    /100;
SCoding_forDNA              = 3*SCoding_forDNA  / (config.sDNwin*2) /100;
SCoding_forDNB              = 3*SCoding_forDNB  / (config.sDNwin/2) /100;
% NS substitutions per cM  (sDNwin, sDNwin*2, sDNwin/2)
SDn                         = SdNcoding .*SCoding_forDN;
SDnA                        = SdNcodingA.*SCoding_forDNA;
SDnB                        = SdNcodingB.*SCoding_forDNB;
% sites in cM  (sDNwin)
[~,~,~,~,~,SSites_forDN   ] = smooth01NonUniform( allsites_gx, ones(size(allsites_gx)), poly_gx, config.sDNwin);
SSites_forDN                = SSites_forDN      / config.sDNwin /100;
% density of coding sites,  sDNwin
SCodingDensity_forDN        = SCoding_forDN    ./ SSites_forDN;
% probability for NS substitution at a SITE,  sDNwin
SdN                         = SdNcoding.*SCodingDensity_forDN;

% % TARGET:   SDn2
% % probability for NS substitution at a CODING SITE  (sDNwin2)
% SdNcoding2                  = smooth01NonUniform( DN_gx,       double(Dn),              poly_gx, config.sDNwin2 ) / 3;
% % coding sites in cM  (sDNwin2)
% [~,~,~,~,~,SCoding_forDN2]  = smooth01NonUniform( codons_gx,   ones(size(codons_gx)),   poly_gx, config.sDNwin2);
% SCoding_forDN2              = 3*SCoding_forDN2  / config.sDNwin2 /100;
% % NS substitutions per cM  (sDNwin2)
% SDn2                        = SdNcoding2 .*SCoding_forDN2;
end

% temp.SCoding = SCoding;
% % temp.SSites = SSites;
% temp.SCodingDensity = SCodingDensity;
% temp.SdNcoding = SdNcoding;
% temp.SdN = SdN;
% temp.SDn = SDn;
% temp.recrate = recrate;
% temp.sdiv = sdiv;
% temp.poly = poly;
% temp.div = div;
% temp.poly_px = poly_px;
% temp.poly_gx = poly_gx;



%% bin variation patterns
% config.minimal_bins_for_corr = 100;

% spatial bins
for k=1:length(config.pSpatBin)
  
  jdx = [1:length(poly_px)]; %find(SCodingDensity>median(SCodingDensity));
  
  anals.binned_spat.edges{k} = [1 : config.pSpatBin(k) : CSChrLen(end)+(length(iC)-1)*config.pChrSpacer+config.pSpatBin(k)]'; %[min(poly_px):config.pSpatBin(k):max(poly_px)+config.pSpatBin(k)]';
  %   anals.binned_spat.edges{k} = [min(poly_gx):config.gSpatBin(k):max(poly_gx)+config.gSpatBin(k)]';
  [anals.binned_spat.gPol{k}, anals.binned_spat.gCodons{k}, gX ] = ...
    poolYbyX(poly_px(jdx), poly(jdx),   -1, 1, 'range', anals.binned_spat.edges{k}, 1); %poly_gx
  anals.binned_spat.gPol{k} = anals.binned_spat.gPol{k} ./anals.binned_spat.gCodons{k};
  % anals.binned_spat.gPol2{k}  = poolYbyX(poly_px(jdx), poly(jdx).^2,  -1, 1, 'range', anals.binned_spat.edges{k}, 1 ) ./anals.binned_spat.gCodons{k}; %poly_gx
  % anals.binned_spat.gDiv{k}   = poolYbyX(poly_px(jdx), div(jdx),      -1, 1, 'range', anals.binned_spat.edges{k}, 1 ) ./anals.binned_spat.gCodons{k}; %poly_gx
  anals.binned_spat.gSDiv{k}  = poolYbyX(poly_px(jdx), sdiv(jdx),     -1, 1, 'range', anals.binned_spat.edges{k}, 1 ) ./anals.binned_spat.gCodons{k}; %poly_gx
  anals.binned_spat.gRat{k}   = poolYbyX(poly_px(jdx), recrate(jdx),  -1, 1, 'range', anals.binned_spat.edges{k}, 1 ) ./anals.binned_spat.gCodons{k}; %poly_gx
  if isfield(fdata{c}, 'Div')
  anals.binned_spat.gSDn{k}   = poolYbyX(poly_px(jdx), SDn(jdx),      -1, 1, 'range', anals.binned_spat.edges{k}, 1 ) ./anals.binned_spat.gCodons{k}; %poly_gx
  anals.binned_spat.gSdN{k}   = poolYbyX(poly_px(jdx), SdN(jdx),      -1, 1, 'range', anals.binned_spat.edges{k}, 1 ) ./anals.binned_spat.gCodons{k}; %poly_gx
  anals.binned_spat.gSCN{k}   = poolYbyX(poly_px(jdx), SCoding(jdx),  -1, 1, 'range', anals.binned_spat.edges{k}, 1 ) ./anals.binned_spat.gCodons{k}; %poly_gx
  anals.binned_spat.gSCD{k}   = poolYbyX(poly_px(jdx), SCodingDensity(jdx),  -1, 1, 'range', anals.binned_spat.edges{k}, 1 ) ./anals.binned_spat.gCodons{k}; %poly_gx
  end

  [anals.binned_spat.stats.C_P2DvC(k),   anals.binned_spat.stats.vol_P2DvC(k),   te,                                 anals.binned_spat.stats.p1_P2DvC(k,:)  ] = corrBinned(anals.binned_spat.gCodons{k}, anals.binned_spat.gRat{k}, anals.binned_spat.gPol{k}./anals.binned_spat.gSDiv{k}, 'Spearman', config.minimal_codons_per_bin);
  [anals.binned_spat.stats.C_P2DvSDn(k), anals.binned_spat.stats.vol_P2DvSDn(k), te,                                 anals.binned_spat.stats.p1_P2DvSDn(k,:)] = corrBinned(anals.binned_spat.gCodons{k}, anals.binned_spat.gSDn{k}, anals.binned_spat.gPol{k}./anals.binned_spat.gSDiv{k}, 'Spearman', config.minimal_codons_per_bin);
  [anals.binned_spat.stats.C_P2DvSCD(k), anals.binned_spat.stats.vol_P2DvSCD(k), te,                                 anals.binned_spat.stats.p1_P2DvSCD(k,:)] = corrBinned(anals.binned_spat.gCodons{k}, anals.binned_spat.gSCD{k}, anals.binned_spat.gPol{k}./anals.binned_spat.gSDiv{k}, 'Spearman', config.minimal_codons_per_bin);
%   [anals.binned_spat.stats.C_P2DvC(k),   anals.binned_spat.stats.vol_P2DvC(k),   anals.binned_spat.stats.qtl_P2DvC(k,:),   anals.binned_spat.stats.p1_P2DvC(k,:)  ] = corrBinned(anals.binned_spat.gCodons{k}, anals.binned_spat.gRat{k}./anals.binned_spat.gCodons{k}, anals.binned_spat.gPol{k}./anals.binned_spat.gSDiv{k}, 'R2');
%   [anals.binned_spat.stats.C_P2DvSDn(k), anals.binned_spat.stats.vol_P2DvSDn(k), anals.binned_spat.stats.qtl_P2DvSDn(k,:), anals.binned_spat.stats.p1_P2DvSDn(k,:)] = corrBinned(anals.binned_spat.gCodons{k}, anals.binned_spat.gSDn{k}./anals.binned_spat.gCodons{k}, anals.binned_spat.gPol{k}./anals.binned_spat.gSDiv{k}, 'R2');
%   [anals.binned_spat.stats.C_P2DvDn(k), anals.binned_spat.stats.vol_P2DvSDn(k)] = corrBinned(anals.binned_spat.gCodons{k}, anals.binned_spat.gPol{k}./anals.binned_spat.gSDiv{k}, anals.binned_spat.gDn{k}./anals.binned_spat.gCodons{k});
  
  idx  = find(anals.binned_spat.gCodons{k}>0);
  % idx2 = find(anals.binned_spat.gDiv{k}>0);
  idx3 = find(anals.binned_spat.gSDiv{k}>0);
  anals.gwa.P(k)   = mean(anals.binned_spat.gPol{k}(idx)); %./anals.binned_spat.gCodons{k}(idx)
  %   anals.gwa.D(k)   = mean(anals.binned_spat.gDiv{k}(idx)); %./anals.binned_spat.gCodons{k}(idx)
  anals.gwa.SD(k)  = mean(anals.binned_spat.gSDiv{k}(idx)); %./anals.binned_spat.gCodons{k}(idx)
  % anals.gwa.P2D(k) = mean(anals.binned_spat.gPol{k}(idx2) ./anals.binned_spat.gDiv{k}(idx2));
  anals.gwa.P2SD(k)= mean(anals.binned_spat.gPol{k}(idx3) ./anals.binned_spat.gSDiv{k}(idx3));
end

anals.gwa.gwhgP   = [0:0.01:0.5];
anals.gwa.gwhgSD  = [0:0.01:0.5];
anals.gwa.gwhgP2SD= [0:0.01:0.5];
anals.gwa.gwhgRec = [0:0.1:5];
anals.gwa.gwhgSDn = [0:10^2:1.5*10^4];
anals.gwa.gwhgSdN = [0:0.001:0.02];
anals.gwa.gwhgSCN = [0:10^4:10^6];
anals.gwa.gwhgSCD = [0:0.1:1];
anals.gwa.gwhP   = histc(poly,            anals.gwa.gwhgP);
anals.gwa.gwhSD  = histc(sdiv,            anals.gwa.gwhgSD);
anals.gwa.gwhP2SD= histc(poly./sdiv,      anals.gwa.gwhgP2SD);
anals.gwa.gwhRec = histc(recrate,         anals.gwa.gwhgRec);
if isfield(fdata{c}, 'Div')
anals.gwa.gwhSDn = histc(SDn,             anals.gwa.gwhgSDn);
anals.gwa.gwhSdN = histc(SdN,             anals.gwa.gwhgSdN);
anals.gwa.gwhSCN = histc(SCoding,         anals.gwa.gwhgSCN);
anals.gwa.gwhSCD = histc(SCodingDensity,  anals.gwa.gwhgSCD);
end


% property bins (RecRate, Dn)
% config.corr_by_prop_bins = [100 25 400 1600];
% config.tail_widths = [0.1 0.05 0.025 0.01 0.005 0.0025];

for b=1:length(config.corr_by_prop_bins)
  ctail_widths{b} = 2*ceil( config.tail_widths*config.corr_by_prop_bins(b)/2 );
end

idx = find( ~isnan(recrate) ); % & SDn>median(SDn)
if isfield(fdata{c}, 'Div')
jdx = [1:length(SDn)]; %find(                   SDn>median(SDn) ); %
end
for b=1:length(config.corr_by_prop_bins)
  
  for p=1:4%5
    
    if p==1
      iidx  = idx;
    else
      iidx  = jdx;
    end
    
    prop  = []; prop2 = [];
    if isfield(fdata{c}, 'Div')
    switch p
      case 1
        prop     = recrate;
        prop2{1} = SCoding;
        prop2{2} = SCodingDensity;
      case 2
        prop     = SDn;
        prop2{1} = SCoding;
        prop2{2} = SCodingDensity;
      case 3
        prop     = SdN;
      case 4
        prop     = SCoding;
    end
    else
    switch p
      case 1
        prop     = recrate;
    end
    end
    
    if ~isempty(prop)
    [binned_prop.gPol, binned_prop.gCodons, ~, binned_prop.edges ] = ...
      poolYbyX(prop(iidx), poly(iidx),   config.corr_by_prop_bins(b), 1, 'percentiles', -1, 1);
    binned_prop.gPol   = binned_prop.gPol ./ binned_prop.gCodons;
    % binned_prop.gPol2  = poolYbyX( prop(iidx), poly(iidx).^2,-1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    % binned_prop.gDiv   = poolYbyX( prop(iidx), div(iidx),    -1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    binned_prop.gSDiv  = poolYbyX( prop(iidx), sdiv(iidx),   -1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    binned_prop.gRat   = poolYbyX( prop(iidx), recrate(iidx),-1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    if ~isempty(prop2)
    binned_prop.gSdN   = poolYbyX( prop(iidx), SdN(iidx),    -1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    binned_prop.gSDn   = poolYbyX( prop(iidx), SDn(iidx),    -1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    binned_prop.gSCN   = poolYbyX( prop(iidx), SCoding(iidx),-1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    binned_prop.gSCD   = poolYbyX( prop(iidx), SCodingDensity(iidx), -1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    end
    %     binned_prop.gProp  = poolYbyX( prop(iidx), prop(iidx),   -1, 1, 'range', binned_prop.edges, 1 ) ./binned_prop.gCodons;
    
    hdx = find(~isnan(binned_prop.gPol)&~isinf(binned_prop.gPol)&~isnan(binned_prop.gSDiv)&~isinf(binned_prop.gSDiv)&(binned_prop.gSDiv~=0));
    midL = floor(length(hdx)/2);
    %     binned_prop.stats.rankcorr = corr( [1:length(hdx)]', binned_prop.gPol(hdx)'./binned_prop.gSDiv(hdx)', 'type', 'Spearman' );
    for t=1:length(config.tail_widths)
      binned_prop.stats.otails(1,t) = mean( binned_prop.gPol(hdx(1:ctail_widths{b}(t))) ./ binned_prop.gSDiv(hdx(1:ctail_widths{b}(t))) );
      binned_prop.stats.otails(2,t) = mean( binned_prop.gPol(hdx(end-ctail_widths{b}(t)+1:end)) ./ binned_prop.gSDiv(hdx(end-ctail_widths{b}(t)+1:end)) );
      binned_prop.stats.otails(3,t) = mean( binned_prop.gPol(hdx) ./ binned_prop.gSDiv(hdx) );
      %     binned_prop.stats.otails(3,t) = mean( binned_prop.gPol(hdx(midL+1-ctail_widths{b}(t)/2:midL+ctail_widths{b}(t)/2)) ./ binned_prop.gSDiv(hdx(midL+1-ctail_widths{b}(t)/2:midL+ctail_widths{b}(t)/2)) );
    end
    
    binned_prop2={[],[]};
    if b==2 & ~isempty(prop2)
      for y=1:length(prop2)
        [binned_prop2{y}.gPol, binned_prop2{y}.gCodons, binned_prop2{y}.edges1, binned_prop2{y}.edges2] = ...
          poolYbyX2D( [prop(iidx) prop2{y}(iidx)], poly(iidx), [config.corr_by_prop_bins(b) 2], 0, {'percentiles', 'percentiles'}, [], 1 );
        binned_prop2{y}.gPol   = binned_prop2{y}.gPol ./binned_prop2{y}.gCodons;
        binned_prop2{y}.gSDiv  = poolYbyX2D( [prop(iidx) prop2{y}(iidx)], sdiv(iidx),    [-1 -1],                    1, {'range',       'range'},       {binned_prop2{y}.edges1, binned_prop2{y}.edges2}, 1 ) ./binned_prop2{y}.gCodons;
        binned_prop2{y}.gRat   = poolYbyX2D( [prop(iidx) prop2{y}(iidx)], recrate(iidx), [-1 -1],                    1, {'range',       'range'},       {binned_prop2{y}.edges1, binned_prop2{y}.edges2}, 1 ) ./binned_prop2{y}.gCodons;
        binned_prop2{y}.gSDn   = poolYbyX2D( [prop(iidx) prop2{y}(iidx)], SDn(iidx),     [-1 -1],                    1, {'range',       'range'},       {binned_prop2{y}.edges1, binned_prop2{y}.edges2}, 1 ) ./binned_prop2{y}.gCodons;
        binned_prop2{y}.gSDnA  = poolYbyX2D( [prop(iidx) prop2{y}(iidx)], SDnA(iidx),    [-1 -1],                    1, {'range',       'range'},       {binned_prop2{y}.edges1, binned_prop2{y}.edges2}, 1 ) ./binned_prop2{y}.gCodons;
        binned_prop2{y}.gSDnB  = poolYbyX2D( [prop(iidx) prop2{y}(iidx)], SDnB(iidx),    [-1 -1],                    1, {'range',       'range'},       {binned_prop2{y}.edges1, binned_prop2{y}.edges2}, 1 ) ./binned_prop2{y}.gCodons;
        binned_prop2{y}.gSdN   = poolYbyX2D( [prop(iidx) prop2{y}(iidx)], SdN(iidx),     [-1 -1],                    1, {'range',       'range'},       {binned_prop2{y}.edges1, binned_prop2{y}.edges2}, 1 ) ./binned_prop2{y}.gCodons;
        binned_prop2{y}.gCN    = poolYbyX2D( [prop(iidx) prop2{y}(iidx)], SCoding(iidx), [-1 -1],                    1, {'range',       'range'},       {binned_prop2{y}.edges1, binned_prop2{y}.edges2}, 1 ) ./binned_prop2{y}.gCodons;
        binned_prop2{y}.gCD    = poolYbyX2D( [prop(iidx) prop2{y}(iidx)], SCodingDensity(iidx), [-1 -1],             1, {'range',       'range'},       {binned_prop2{y}.edges1, binned_prop2{y}.edges2}, 1 ) ./binned_prop2{y}.gCodons;
      end
    end
    end
    
    if isfield(fdata{c}, 'Div')
    switch p
      case 1
        anals.binned_prop.c{b}   = binned_prop;
        anals.binned_prop.cCN{b} = binned_prop2{1};
        anals.binned_prop.cCD{b} = binned_prop2{2};
      case 2
        anals.binned_prop.SDn{b}   = binned_prop;
        anals.binned_prop.SDnCN{b} = binned_prop2{1};
        anals.binned_prop.SDnCD{b} = binned_prop2{2};
      case 3
        anals.binned_prop.SdN{b}   = binned_prop;
      case 4
        anals.binned_prop.SCN{b}   = binned_prop;
    end
    else
    switch p
      case 1
        anals.binned_prop.c{b}   = binned_prop;
    end
    end
  end
  
  
  
  
  % prepare for fit to the functional form Wiehe&Kim&Innan&Stephan
  if isfield(fdata{c}, 'Div')
  dd_byC{b}   = anals.binned_prop.c{  b}.gPol ./anals.binned_prop.c{  b}.gSDiv;
  dd_bySDn{b} = anals.binned_prop.SDn{b}.gPol ./anals.binned_prop.SDn{b}.gSDiv;
  cc{b}       = anals.binned_prop.c{  b}.gRat * 10^-8;
  dN_byC{b}   = anals.binned_prop.c{b}.gSdN;
  dN{b}       = anals.binned_prop.SDn{b}.gSdN;
  bDn{b}      = anals.binned_prop.SDn{b}.gSDn;
  cc_bySDn{b} = anals.binned_prop.SDn{b}.gRat * 10^-8;
  jjj{b}      = ~isnan(cc{b});
  iii{b}      = ~isnan(dN{b});
  EdN{b}      = ones(size(dd_byC{  b})) * sum(dN{b}(iii{b}).*anals.binned_prop.SDn{b}.gCodons(iii{b})) ./ sum(anals.binned_prop.SDn{b}.gCodons(iii{b}));
  Ecc{b}      = ones(size(dd_bySDn{b})) * sum(cc{b}(jjj{b}).*anals.binned_prop.c{  b}.gCodons(jjj{b})) ./ sum(anals.binned_prop.c{  b}.gCodons(jjj{b}));
  end
end


if isfield(fdata{c}, 'Div')
  
b = 1;
if isempty(toymodelfits) 
  toymodelfits.fitBGS_byC         = fitBGSnRHHbyCnDn(  cc{b},              dN_byC{b},  dd_byC{  b}, 'BGS', jjj{b} ); %dN_byC{b} - dN doesn't play any role for this model!
  toymodelfits.fitRHH_byC         = fitBGSnRHHbyCnDn(  cc{b},              dN_byC{b},  dd_byC{  b}, 'RHH', jjj{b} );
  %  toymodelfits.fitRHH_byC         = fitBGSnRHHbyCnDn(  cc{b},              EdN{b},  dd_byC{  b}, 'RHH', jjj{b} );
  [toymodelfits.fitRHH_bySDn, idx1]= fitBGSnRHHbyCnDn(  cc_bySDn{b},        dN{b},   dd_bySDn{b}, 'RHH2', iii{b} );
  [toymodelfits.fitRHH_bySDn2,idx2]= fitBGSnRHHbyCnDn(  ones(size(bDn{b})), bDn{b},  dd_bySDn{b}, 'RHH2', iii{b} );
  toymodelfits.fitRHH_bySDn.SDn   = anals.binned_prop.SDn{b}.gSDn(idx1);
  toymodelfits.fitRHH_bySDn2.SDn  = anals.binned_prop.SDn{b}.gSDn(idx2);
  %   toymodelfits.fitRHH_bySDn = fitBGSnRHHbyCnDn(  cc_bySDn{b},  dN{b}, dd_bySDn{b}, 'RHH', iii{b} );
end
anals.toymodelfits = toymodelfits;
%
for b=1:length(config.corr_by_prop_bins)
  [te, anals.binned_prop.c{b}.fitBGS_byC]      = serrBGSnRHH(anals.toymodelfits.fitBGS_byC.params,     cc{b},       dN_byC{b},  dd_byC{b}, jjj{b});
  % [te, anals.binned_prop.c{b}.fitBGS_byC]      = serrBGSnRHH(anals.toymodelfits.fitBGS_byC.params,     cc{b},       EdN{b},     dd_byC{b}, jjj{b});
  [te, anals.binned_prop.c{b}.fitRHH_byC]      = serrBGSnRHH(anals.toymodelfits.fitRHH_byC.params,     cc{b},       dN_byC{b},  dd_byC{b}, jjj{b});
  % [te, anals.binned_prop.c{b}.fitRHH_byC]      = serrBGSnRHH(anals.toymodelfits.fitRHH_byC.params,     cc{b},       EdN{b},     dd_byC{b}, jjj{b});
  [te, anals.binned_prop.SDn{b}.fitRHH_bySDn]  = serrBGSnRHH(anals.toymodelfits.fitRHH_bySDn.params,   cc_bySDn{b}, dN{b},      dd_bySDn{b}, iii{b});
  % [te, anals.binned_prop.SDn{b}.fitRHH_bySDn]  = serrBGSnRHH(anals.toymodelfits.fitRHH_bySDn.params,   ones(size(cc{b})), bDn{b}, dd_bySDn{b}, iii{b});
  anals.binned_prop.c{b}.stats.R2_BGS_byC     = Rsquare( dd_byC{b}( jjj{b}),  anals.binned_prop.c{b}.fitBGS_byC.pred( jjj{b}) );
  anals.binned_prop.c{b}.stats.R2_RHH_byC     = Rsquare( dd_byC{b}( jjj{b}),  anals.binned_prop.c{b}.fitRHH_byC.pred( jjj{b}) );
  anals.binned_prop.SDn{b}.stats.R2_RHH_bySDn = Rsquare( dd_bySDn{b}(iii{b}), anals.binned_prop.SDn{b}.fitRHH_bySDn.pred( iii{b}) );
  %   anals.binned_prop.SDn{b}.stats.R2_RHH_bySDn = Rsquare( dd_bySDn{b}(iii{b}), anals.binned_prop.SDn{b}.fitRHH_bySDn.pred );
  
  anals.binned_prop.SDn{b}.stats.RedSW_RHH_bySDn = mean(anals.binned_prop.SDn{b}.fitRHH_bySDn.pred(iii{b}))/anals.toymodelfits.fitRHH_bySDn.params(1);
  anals.binned_prop.c{b}.stats.RedSW_RHH_byC     = mean(anals.binned_prop.c{b  }.fitRHH_byC.pred(jjj{b})  )/anals.toymodelfits.fitRHH_byC.params(1);
  anals.binned_prop.c{b}.stats.RedBS_BGS_byC     = mean(anals.binned_prop.c{b  }.fitBGS_byC.pred(jjj{b})  )/anals.toymodelfits.fitBGS_byC.params(1);
  
  idx = ~isnan(anals.binned_prop.c{b}.fitBGS_byC.pred) & ~isnan(anals.binned_prop.c{b}.fitBGS_byC.obs);
  anals.binned_prop.c{b}.stats.rankcorr    = corr( anals.binned_prop.c{b}.fitBGS_byC.pred(idx)', anals.binned_prop.c{b}.fitBGS_byC.obs(idx)', 'type', 'Spearman' );
  anals.binned_prop.c{b}.stats.lincorr     = corr( anals.binned_prop.c{b}.fitBGS_byC.pred(idx)', anals.binned_prop.c{b}.fitBGS_byC.obs(idx)', 'type', 'Pearson' );
  idx = ~isnan(anals.binned_prop.c{b}.fitRHH_byC.pred) & ~isnan(anals.binned_prop.c{b}.fitRHH_byC.obs);
  anals.binned_prop.c{b}.stats.rankcorr2   = corr( anals.binned_prop.c{b}.fitRHH_byC.pred(idx)', anals.binned_prop.c{b}.fitRHH_byC.obs(idx)', 'type', 'Spearman' );
  anals.binned_prop.c{b}.stats.lincorr2    = corr( anals.binned_prop.c{b}.fitRHH_byC.pred(idx)', anals.binned_prop.c{b}.fitRHH_byC.obs(idx)', 'type', 'Pearson' );
  idx = ~isnan(anals.binned_prop.SDn{b}.fitRHH_bySDn.pred) & ~isnan(anals.binned_prop.SDn{b}.fitRHH_bySDn.obs);
  anals.binned_prop.SDn{b}.stats.rankcorr  = corr( anals.binned_prop.SDn{b}.fitRHH_bySDn.pred(idx)', anals.binned_prop.SDn{b}.fitRHH_bySDn.obs(idx)', 'type', 'Spearman' );
  anals.binned_prop.SDn{b}.stats.lincorr   = corr( anals.binned_prop.SDn{b}.fitRHH_bySDn.pred(idx)', anals.binned_prop.SDn{b}.fitRHH_bySDn.obs(idx)', 'type', 'Pearson' );
  
  binned_prop = anals.binned_prop.c{b};
  hdx = find(~isnan(binned_prop.gPol)&~isinf(binned_prop.gPol)&~isnan(binned_prop.gSDiv)&~isinf(binned_prop.gSDiv)&(binned_prop.gSDiv~=0));
  midL = floor(length(hdx)/2);
  for t=1:length(config.tail_widths)
    binned_prop.stats.ptails(1,t) = mean( anals.binned_prop.c{b}.fitBGS_byC.pred(hdx(1:ctail_widths{b}(t))) );
    binned_prop.stats.ptails(2,t) = mean( anals.binned_prop.c{b}.fitBGS_byC.pred(hdx(end-ctail_widths{b}(t)+1:end)) );
    binned_prop.stats.ptails(3,t) = mean( anals.binned_prop.c{b}.fitBGS_byC.pred(hdx) );
  end
  anals.binned_prop.c{b}.stats.ptails_BGS_byC = binned_prop.stats.ptails;
  
  binned_prop = anals.binned_prop.c{b};
  hdx = find(~isnan(binned_prop.gPol)&~isinf(binned_prop.gPol)&~isnan(binned_prop.gSDiv)&~isinf(binned_prop.gSDiv)&(binned_prop.gSDiv~=0));
  midL = floor(length(hdx)/2);
  for t=1:length(config.tail_widths)
    binned_prop.stats.ptails(1,t) = mean( anals.binned_prop.c{b}.fitRHH_byC.pred(hdx(1:ctail_widths{b}(t))) );
    binned_prop.stats.ptails(2,t) = mean( anals.binned_prop.c{b}.fitRHH_byC.pred(hdx(end-ctail_widths{b}(t)+1:end)) );
    binned_prop.stats.ptails(3,t) = mean( anals.binned_prop.c{b}.fitRHH_byC.pred(hdx) );
  end
  anals.binned_prop.c{b}.stats.ptails_RHH_byC = binned_prop.stats.ptails;
  
  binned_prop = anals.binned_prop.SDn{b};
  hdx = find(~isnan(binned_prop.gPol)&~isinf(binned_prop.gPol)&~isnan(binned_prop.gSDiv)&~isinf(binned_prop.gSDiv)&(binned_prop.gSDiv~=0));
  midL = floor(length(hdx)/2);
  for t=1:length(config.tail_widths)
    binned_prop.stats.ptails(1,t) = mean( anals.binned_prop.SDn{b}.fitRHH_bySDn.pred(hdx(1:ctail_widths{b}(t))) );
    binned_prop.stats.ptails(2,t) = mean( anals.binned_prop.SDn{b}.fitRHH_bySDn.pred(hdx(end-ctail_widths{b}(t)+1:end)) );
    binned_prop.stats.ptails(3,t) = mean( anals.binned_prop.SDn{b}.fitRHH_bySDn.pred(hdx) );
  end
  anals.binned_prop.SDn{b}.stats.ptails_RHH_bySDn = binned_prop.stats.ptails;
end

end

% these parameters are enforced to ensure a reasonable collation processing time
% config.gdelta_COL = 10^-5; %10^-6 10^-8
% config.halfway = 0; %100;
% config.gL_COL = 0.11*0.01; %0.5*0.01;
% config.sort_strand_dir = 1;
% config.COL_R2_window = floor(min( [1 2 0.5]*110, 1+2*config.gL_COL/config.gdelta_COL ));

% PATCH
for c=iC
  config.pos{c} = fdata{c}.Poly.pos(preCalc{1}{c}.idxc);
end

%% collated observed diversity and divergence
if calculate>1
  
  if ~isfield(config.collate, 'focals_idx') || isempty(config.collate.focals_idx)
    for c=iC
      config.collate.focals_idx.NS{ c} = [];
      config.collate.focals_idx.SYN{c} = [];
      config.collate.focals_idx.UTR{c} = [];
    end
  end
  
  anals.collated{1}.NS      = collateCodons( fdata, focals.SubsNS,   config.collate, [], config.collate.focals_idx.NS,  config.predef_idx_test ); %config.Col
  anals.collated{1}.SYN     = collateCodons( fdata, focals.SubsSYN,  config.collate, [], config.collate.focals_idx.SYN, config.predef_idx_test ); %config.Col
  anals.collated{1}.UTR     = collateCodons( fdata, focals.SubsUTR,  config.collate, [], config.collate.focals_idx.UTR, config.predef_idx_test ); %config.Col
  %   anals.collated{1}.INT     = collateCodons( fdata(iC), focals.SubsInt(iC),  config ); %config.Col
  %   anals.collated{1}.Exons   = collateCodons( fdata(iC), focals.ExonTips(iC), config ); %config.Col
  
  configZ = config.collate;
  configZ.gL_COL     = configZ.gL_COL/50;
  configZ.gdelta_COL = configZ.gdelta_COL/100;
  anals.collated{2}.NS      = collateCodons( fdata, focals.SubsNS,   configZ, [], config.collate.focals_idx.NS,  config.predef_idx_test ); %config.Col
  anals.collated{2}.SYN     = collateCodons( fdata, focals.SubsSYN,  configZ, [], config.collate.focals_idx.SYN, config.predef_idx_test ); %config.Col
  anals.collated{2}.UTR     = collateCodons( fdata, focals.SubsUTR,  configZ, [], config.collate.focals_idx.UTR, config.predef_idx_test ); %config.Col
  %   anals.collated{2}.INT     = collateCodons( fdata(iC), focals.SubsInt(iC),  configZ ); %config.Col
  %   anals.collated{2}.Exons   = collateCodons( fdata(iC), focals.ExonTips(iC), configZ ); %config.Col
  
  

  if isfield(fdata{c}, 'Div')
  % quantify clustering of subs, codons
  for t=1:2
    switch t
      case 1
        cur_c_coll = anals.collated{1}.NS;
      case 2
        cur_c_coll = anals.collated{1}.SYN;
      case 3
        cur_c_coll = anals.collated{1}.UTR;
    end
    if isfield(cur_c_coll, 'clDn')
      
      LL = (length(cur_c_coll.Hs.radius)-1)/2;
      
      % subs per bin, corrected for missing data
      ww1 = cur_c_coll.clDn.property / cur_c_coll.clDn.focals;
      % subs per cM, corrected for missing data
      ww  = ww1 / (100*cur_c_coll.Hs.gbin);
      % the radius of the excess peak is defined as the distance between points where the level is 110% of the median value which represents the baselevel
      if t==1
        excess_subs_med    = prctile(ww, 50);
        excess_subs_th     = 1.1*excess_subs_med;
        excess_subs_radius = find( 0.5*(sort(ww(1:LL),'descend')+sort(ww(LL+2:end),'descend')) < excess_subs_th, 1 ) - 1;
      end
      cluststats.excess_subs_med    = excess_subs_med;
      cluststats.excess_subs_th     = excess_subs_th;
      cluststats.excess_subs_radius = excess_subs_radius;
      cluststats.excess_subs_dist = 100*cur_c_coll.Hs.gbin * cluststats.excess_subs_radius;
      % the excess of subs within the radius is the sum of their excess above the median level
      cluststats.excess_subs = sum( ww1(1+LL-cluststats.excess_subs_radius:1+LL+cluststats.excess_subs_radius) - excess_subs_med );
      
      % codons per bin
      vv1 = cur_c_coll.clDnAnno.codons / cur_c_coll.clDnAnno.focals;
      % codons per cM
      vv  = vv1 / (100*cur_c_coll.Hs.gbin);
      % the radius of the excess peak is defined as the distance between points where the level is 110% of the mdeian value which represents the baselevel
      if t==1
        excess_codons_med    = prctile(vv, 50);
        excess_codons_th     = 1.1*excess_codons_med;
        excess_codons_radius = find( 0.5*(sort(vv(1:LL),'descend')+sort(vv(LL+2:end),'descend')) < excess_codons_th, 1 ) - 1;
      end
      cluststats.excess_codons_med    = excess_codons_med;
      cluststats.excess_codons_th     = excess_codons_th;
      cluststats.excess_codons_radius = excess_codons_radius;
      cluststats.excess_codons_dist = 100*cur_c_coll.Hs.gbin * cluststats.excess_codons_radius;
      % the excess of codons within the radius is the sum of their excess above the median level
      cluststats.excess_codons = sum( vv1(1+LL-cluststats.excess_codons_radius:1+LL+cluststats.excess_codons_radius) - excess_codons_med );
      
      switch t
        case 1
          anals.collated{1}.NS.stats  = cluststats;
        case 2
          anals.collated{1}.SYN.stats = cluststats;
        case 3
          anals.collated{1}.UTR.stats = cluststats;
      end
    end
  end
  end
end


anals.config = config;




end

