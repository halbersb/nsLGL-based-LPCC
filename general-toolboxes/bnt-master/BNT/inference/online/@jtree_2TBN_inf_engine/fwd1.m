function [f, logscale] = fwd1(engine, ev, t)
% Forwards pass for slice 1.

bnet = bnet_from_engine(engine);
ss = bnet.nnodes_per_slice;

CPDpot = cell(1,ss);
for n=1:ss
  fam = family(bnet.dag, n, 1);
  if size(bnet.equiv_class,2)<size(bnet.equiv_class,1) % dan added the condition
    e = bnet.equiv_class(n, 1);
  else
    e = bnet.equiv_class(1, n);
  end
  CPDpot{n} = convert_to_pot(bnet.CPD{e}, engine.pot_type, fam(:), ev);
end       
f.t = t;
f.evidence = ev;

pots = CPDpot;
slice1 = 1:ss;
CPDclqs = engine.clq_ass_to_node1(slice1);

[f.clpot, f.seppot] =  init_pot(engine.jtree_engine1, CPDclqs, CPDpot, engine.pot_type, engine.observed1);
[f.clpot, f.seppot] = collect_evidence(engine.jtree_engine1, f.clpot, f.seppot);
for c=1:length(f.clpot)
  [f.clpot{c}, ll(c)] = normalize_pot(f.clpot{c});
end
logscale = ll(engine.root1);

