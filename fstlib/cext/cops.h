#ifndef FST_LIBRARY_H_
#define FST_LIBRARY_H_

#include <fst/fstlib.h>
#include <fst/script/fstscript.h>
#include <fst/script/script-impl.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <map>

//DEFINE_int32(v, 10, "v");
const float mydelta = 1.0F/(8192.0F*4);

using namespace fst;

//void shortest_path(script::FstClass &model, script::FstClass &input, script::FstClass &output, script::MutableFstClass* path, std::vector<script::WeightClass>* distance) {
void shortest_path(script::FstClass &model, script::FstClass &input, script::FstClass &output, script::MutableFstClass* path) {	
	const Fst<StdArc> *tfst = model.GetFst<StdArc>();
	const Fst<StdArc> *ifst1 = input.GetFst<StdArc>();
	const Fst<StdArc> *ifst2 = output.GetFst<StdArc>();
	MutableFst<StdArc> *ofst = path->GetMutableFst<StdArc>();

	ArcSortFst<StdArc, OLabelCompare<StdArc> > input_sorted = ArcSortFst<StdArc, OLabelCompare<StdArc> >(*ifst1, OLabelCompare<StdArc>());
	ArcSortFst<StdArc, ILabelCompare<StdArc> > output_sorted = ArcSortFst<StdArc, ILabelCompare<StdArc> >(*ifst2, ILabelCompare<StdArc>());

	// set compose options
	ComposeFstOptions<StdArc> co;
	co.gc_limit=0;

	// Container for composition result.
	ComposeFst<StdArc> middle = ComposeFst<StdArc>(input_sorted, *tfst, co);
	ComposeFst<StdArc> result = ComposeFst<StdArc>(middle, output_sorted, co);

	std::vector<StdArc::Weight> typed_distance;
	//script::WeightClass distance;

	NaturalShortestFirstQueue<StdArc::StateId, StdArc::Weight> state_queue(typed_distance);
	ShortestPathOptions<StdArc, NaturalShortestFirstQueue<StdArc::StateId, StdArc::Weight>, AnyArcFilter<StdArc> > opts(&state_queue, AnyArcFilter<StdArc>(), 1, false, false, mydelta);
	opts.first_path=true;

	ShortestPath(result, ofst, &typed_distance, opts);
	//StdArc::StateId start = path->Start();

	//distance = script::WeightClass(typed_distance[start]);
	//return distance;
}

script::WeightClass shortest_distance_std(script::FstClass &model, script::FstClass &input, script::FstClass &output) {	
	const Fst<StdArc> *tfst = model.GetFst<StdArc>();
	const Fst<StdArc> *ifst1 = input.GetFst<StdArc>();
	const Fst<StdArc> *ifst2 = output.GetFst<StdArc>();
	
	ArcSortFst<StdArc, OLabelCompare<StdArc> > input_sorted = ArcSortFst<StdArc, OLabelCompare<StdArc> >(*ifst1, OLabelCompare<StdArc>());
	ArcSortFst<StdArc, ILabelCompare<StdArc> > output_sorted = ArcSortFst<StdArc, ILabelCompare<StdArc> >(*ifst2, ILabelCompare<StdArc>());

	// set compose options
	ComposeFstOptions<StdArc> co;
	co.gc_limit=0;

	// Container for composition result.
	ComposeFst<StdArc> middle = ComposeFst<StdArc>(input_sorted, *tfst, co);
	ComposeFst<StdArc> result = ComposeFst<StdArc>(middle, output_sorted, co);

	std::vector<StdArc::Weight> typed_distance;
	StdArc::Weight retval;
	script::WeightClass distance;

	NaturalShortestFirstQueue<StdArc::StateId, StdArc::Weight> state_queue(typed_distance);
	ShortestDistanceOptions<StdArc, NaturalShortestFirstQueue<StdArc::StateId, StdArc::Weight>, AnyArcFilter<StdArc> > opts(&state_queue, AnyArcFilter<StdArc>());
	opts.first_path=true;
	
	using StateId = typename StdArc::StateId;
	using Weight = typename StdArc::Weight;

	if (Weight::Properties() & kRightSemiring) {
		ShortestDistance(result, &typed_distance, opts);
		if (typed_distance.size() == 1 && !typed_distance[0].Member()) {
			retval =  StdArc::Weight::NoWeight();
		}
		Adder<Weight> adder;  // maintains cumulative sum accurately
		for (StateId state = 0; state < typed_distance.size(); ++state) {
			adder.Add(Times(typed_distance[state], result.Final(state)));
		}
		retval = adder.Sum();
	} else {
		ShortestDistance(result, &typed_distance, true, mydelta);
		const auto state = result.Start();
		if (typed_distance.size() == 1 && !typed_distance[0].Member()) {
			retval = StdArc::Weight::NoWeight();
		}
		retval = state != kNoStateId && state < typed_distance.size() ? typed_distance[state] : Weight::Zero();
	}
	
	distance = script::WeightClass(retval);
	return distance;
}

script::WeightClass shortest_distance_log(script::FstClass &model, script::FstClass &input, script::FstClass &output) {	
	const Fst<LogArc> *tfst = model.GetFst<LogArc>();
	const Fst<LogArc> *ifst1 = input.GetFst<LogArc>();
	const Fst<LogArc> *ifst2 = output.GetFst<LogArc>();
	
	ArcSortFst<LogArc, OLabelCompare<LogArc> > input_sorted = ArcSortFst<LogArc, OLabelCompare<LogArc> >(*ifst1, OLabelCompare<LogArc>());
	ArcSortFst<LogArc, ILabelCompare<LogArc> > output_sorted = ArcSortFst<LogArc, ILabelCompare<LogArc> >(*ifst2, ILabelCompare<LogArc>());

	// set compose options
	ComposeFstOptions<LogArc> co;
	//co.gc_limit=0;

	// Container for composition result.
	ComposeFst<LogArc> middle = ComposeFst<LogArc>(input_sorted, *tfst, co);
	ComposeFst<LogArc> result = ComposeFst<LogArc>(middle, output_sorted, co);

	std::vector<LogArc::Weight> typed_distance;
	LogArc::Weight retval;
	script::WeightClass distance;

	//NaturalShortestFirstQueue<LogArc::StateId, LogArc::Weight> state_queue(typed_distance);
	//ShortestDistanceOptions<LogArc, NaturalShortestFirstQueue<LogArc::StateId, LogArc::Weight>, AnyArcFilter<LogArc> > opts(&state_queue, AnyArcFilter<LogArc>());
	//opts.first_path=true;
	
	using StateId = typename LogArc::StateId;
	using Weight = typename LogArc::Weight;

	if (Weight::Properties() & kRightSemiring) {
		ShortestDistance(result, &typed_distance);
		if (typed_distance.size() == 1 && !typed_distance[0].Member()) {
			retval =  LogArc::Weight::NoWeight();
		}
		Adder<Weight> adder;  // maintains cumulative sum accurately
		for (StateId state = 0; state < typed_distance.size(); ++state) {
			adder.Add(Times(typed_distance[state], result.Final(state)));
		}
		retval = adder.Sum();
	} else {
		ShortestDistance(result, &typed_distance, true, mydelta);
		const auto state = result.Start();
		if (typed_distance.size() == 1 && !typed_distance[0].Member()) {
			retval = LogArc::Weight::NoWeight();
		}
		retval = state != kNoStateId && state < typed_distance.size() ? typed_distance[state] : Weight::Zero();
	}
	
	distance = script::WeightClass(retval);
	return distance;
}

script::WeightClass kernel_score_std_impl(script::FstClass &model, script::FstClass &input, script::FstClass &output) {	
	const Fst<StdArc> *tfst = model.GetFst<StdArc>();
	const Fst<StdArc> *ifst1 = input.GetFst<StdArc>();
	const Fst<StdArc> *ifst2 = output.GetFst<StdArc>();
	
	ArcSortFst<StdArc, ILabelCompare<StdArc> > input_sorted = ArcSortFst<StdArc, ILabelCompare<StdArc> >(*ifst1, ILabelCompare<StdArc>());
	ArcSortFst<StdArc, ILabelCompare<StdArc> > output_sorted = ArcSortFst<StdArc, ILabelCompare<StdArc> >(*ifst2, ILabelCompare<StdArc>());

	// set compose options
	ComposeFstOptions<StdArc> co;
	//co.gc_limit=0;

	// Container for composition result.
	InvertFst<StdArc> middle1 = InvertFst<StdArc>(ComposeFst<StdArc>(*tfst, input_sorted, co));
	ComposeFst<StdArc> middle2 = ComposeFst<StdArc>(*tfst, output_sorted, co);
	ArcSortFst<StdArc, ILabelCompare<StdArc> > middle2_sorted = ArcSortFst<StdArc, ILabelCompare<StdArc> >(middle2, ILabelCompare<StdArc>());
	//ArcSortFst<StdArc, OLabelCompare<StdArc> > middle1_sorted = ArcSortFst<StdArc, OLabelCompare<StdArc> >(middle1, OLabelCompare<StdArc>());
	ComposeFst<StdArc> result = ComposeFst<StdArc>(middle1, middle2_sorted, co);
	

	std::vector<StdArc::Weight> typed_distance;
	StdArc::Weight retval;
	script::WeightClass distance;

	NaturalShortestFirstQueue<StdArc::StateId, StdArc::Weight> state_queue(typed_distance);
	ShortestDistanceOptions<StdArc, NaturalShortestFirstQueue<StdArc::StateId, StdArc::Weight>, AnyArcFilter<StdArc> > opts(&state_queue, AnyArcFilter<StdArc>());
	opts.first_path=true;
	
	using StateId = typename StdArc::StateId;
	using Weight = typename StdArc::Weight;

	if (Weight::Properties() & kRightSemiring) {
		ShortestDistance(result, &typed_distance, opts);
		if (typed_distance.size() == 1 && !typed_distance[0].Member()) {
			retval =  StdArc::Weight::NoWeight();
		}
		Adder<Weight> adder;  // maintains cumulative sum accurately
		for (StateId state = 0; state < typed_distance.size(); ++state) {
			adder.Add(Times(typed_distance[state], result.Final(state)));
		}
		retval = adder.Sum();
	} else {
		ShortestDistance(result, &typed_distance, true, mydelta);
		const auto state = result.Start();
		if (typed_distance.size() == 1 && !typed_distance[0].Member()) {
			retval = StdArc::Weight::NoWeight();
		}
		retval = state != kNoStateId && state < typed_distance.size() ? typed_distance[state] : Weight::Zero();
	}
	
	distance = script::WeightClass(retval);
	return distance;
}

#endif