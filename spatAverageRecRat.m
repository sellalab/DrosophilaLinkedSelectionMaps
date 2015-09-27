


function grat = spatAverageRecRat( genmap, pos, rec_spat_window )

gpos0 = applyGenmap2pos( genmap, max(min(genmap.pos), pos) );

gpos1 = applyGenmap2pos( genmap, max(min(genmap.pos), pos-rec_spat_window/2) ); %
gpos2 = applyGenmap2pos( genmap, min(max(genmap.pos), pos+rec_spat_window/2) ); %pos+rec_spat_window/2

window = double(max( 1, min(max(genmap.pos), pos+rec_spat_window/2) - max(min(genmap.pos), pos-rec_spat_window/2) ));

grat = 10^8 * (gpos2-gpos1) ./ window;


end

