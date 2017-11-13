/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_VARIANT_BUFFER_HH
#define VCF2MULTIALIGN_VARIANT_BUFFER_HH

#include <boost/container/set.hpp> // For an extract-capable multiset.
#include <mutex>
#include <string>
#include <vcf2multialign/dispatch_fn.hh>
#include <vcf2multialign/types.hh>
#include <vcf2multialign/vcf_reader.hh>


namespace vcf2multialign {
	
	
	struct variant_buffer_delegate
	{
		virtual ~variant_buffer_delegate() {}
		virtual void handle_variant(variant &variant) = 0;
		virtual void finish() = 0;
	};
	
	
	class variant_buffer
	{
		friend void swap(variant_buffer &lhs, variant_buffer &rhs);
		
	protected:
		struct cmp_variant
		{
			bool operator()(variant const &lhs, variant const &rhs) const
			{
				return lhs.ref().size() > rhs.ref().size();
			}
		};
		
		typedef boost::container::multiset <variant, cmp_variant>	variant_set;
		typedef std::vector <variant_set::node_type>				variant_vector;
		
		struct nonmovable
		{
			std::mutex		buffer_mutex{};
			variant_vector	buffer;			// Actually movable but not copyable.
			
			nonmovable() = default;
			nonmovable(nonmovable const &) {}
			nonmovable(nonmovable &&) {}
			nonmovable &operator=(nonmovable const &other) & { return *this; }
			nonmovable &operator=(nonmovable &&other) & { return *this; }
		};
		
	protected:
		nonmovable							m_nm{};
		vcf_reader							*m_reader{};
		variant_buffer_delegate				*m_delegate{};
		dispatch_ptr <dispatch_queue_t>		m_worker_queue{};
		dispatch_ptr <dispatch_semaphore_t>	m_process_sema{};
		variant_set							m_factory;
		variant_set							m_prepared_variants;
		std::size_t							m_previous_pos{};
		
	protected:
		void return_node_to_buffer(variant_set::node_type &&node);
		bool get_node_from_buffer(variant_set::node_type &node);

	public:
		variant_buffer() = default;
		
		variant_buffer(
			vcf_reader &reader,
			dispatch_ptr <dispatch_queue_t> const &worker_queue,
			variant_buffer_delegate &delegate
		):
			m_reader(&reader),
			m_delegate(&delegate),
			m_worker_queue(worker_queue),
			m_process_sema(dispatch_semaphore_create(10), false)
		{
		}
		
		vcf_reader &reader() { return *m_reader; }
		void read_input();
		void process_input(variant_set &variants);
		void set_vcf_reader(vcf_reader &reader) { m_reader = &reader; }
		void set_delegate(variant_buffer_delegate &delegate) { m_delegate = &delegate; }
	};
}

#endif
