/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <boost/range/adaptor/map.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/function_output_iterator.hpp>
#include <vcf2multialign/reduced_subgraph.hh>


namespace adapt = boost::adaptors;
namespace v2m = vcf2multialign;


namespace {
	
	struct tagged_sample_id
	{
		v2m::sample_id						sample_id;
		v2m::reduced_subgraph const			*target_subgraph{nullptr};
		v2m::reduced_subgraph::path_index	target_path_index{};
		
		tagged_sample_id(
			v2m::sample_id const &sample_id_,
			v2m::reduced_subgraph const &target_subgraph_,
			v2m::reduced_subgraph::path_index const target_path_index_
		):
			sample_id(sample_id_),
			target_subgraph(&target_subgraph_),
			target_path_index(target_path_index_)
		{
		}
		
		bool operator<(tagged_sample_id const &rhs) const
		{
			return sample_id < rhs.sample_id;
		}
	};
	
	
	class tag_sample_id
	{
	public:
		typedef tagged_sample_id result_type;
		
	protected:
		v2m::reduced_subgraph const			*m_target_subgraph{nullptr};
		v2m::reduced_subgraph::path_index	m_target_path_index{};
		
	public:
		tag_sample_id() = default;
		
		tag_sample_id(
			v2m::reduced_subgraph const &target_subgraph,
			v2m::reduced_subgraph::path_index const target_path_index
		):
			m_target_subgraph(&target_subgraph),
			m_target_path_index(target_path_index)
		{
		}
		
		tagged_sample_id operator()(v2m::sample_id const &sample_id) const
		{
			return tagged_sample_id(sample_id, *m_target_subgraph, m_target_path_index);
		}
	};
	
	
	class check_sample_id
	{
	protected:
		std::size_t *m_weight{nullptr};
		
	public:
		check_sample_id() = default;
		
		check_sample_id(std::size_t &weight): m_weight(&weight) {}
		
		void operator()(tagged_sample_id const &tagged_sample_id) const
		{
			auto const &sample_id(tagged_sample_id.sample_id);
			auto const &graph(*tagged_sample_id.target_subgraph);
			auto const path_idx(tagged_sample_id.target_path_index);
			
			v2m::reduced_subgraph::sequence_index sample_seq_idx{};
			auto const st_1(graph.sample_sequence_index(sample_id, sample_seq_idx));
			v2m::always_assert(st_1);
			
			v2m::reduced_subgraph::sequence_index path_seq_idx{};
			auto const st_2(graph.path_sequence_index(path_idx, path_seq_idx));
			v2m::always_assert(st_2);
			
			if (sample_seq_idx != path_seq_idx)
				++*m_weight;
		}
	};
}


namespace vcf2multialign {
	
	std::size_t edge_weight(
		reduced_subgraph const &lhs,
		reduced_subgraph const &rhs,
		reduced_subgraph::path_index const li,
		reduced_subgraph::path_index const ri
	)
	{
		// std::set_symmetric_difference does not mark the origin of each element written to the output iterator,
		// while we would like to check the equivalence classes of those sample ids w.r.t. to their ALT index sequences.
		// We use check_sample_id with a function output iterator to determine and count those but before doing that,
		// we convert the iterator pairs into ranges and pipe them to range adaptors that tag the sample ids with
		// the subgraphs.
		auto const lhs_samples_it_pair(lhs.path_samples(li));
		auto const rhs_samples_it_pair(rhs.path_samples(ri));
		auto const lhs_range(boost::make_iterator_range(lhs_samples_it_pair.first, lhs_samples_it_pair.second));
		auto const rhs_range(boost::make_iterator_range(rhs_samples_it_pair.first, rhs_samples_it_pair.second));
		
		// Tag with target subgraph and path.
		tag_sample_id tag_lhs(rhs, ri);
		tag_sample_id tag_rhs(lhs, li);
		std::size_t weight(0);
		
		boost::set_symmetric_difference(
			lhs_range | adapt::map_values | adapt::transformed(tag_lhs),
			rhs_range | adapt::map_values | adapt::transformed(tag_rhs),
			boost::make_function_output_iterator(check_sample_id(weight))
		);
		
		return weight;
	}
}
