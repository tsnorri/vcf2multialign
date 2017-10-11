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
		
		// Movable instance variables.
		struct data
		{
			vcf_reader							*m_reader{};
			variant_buffer_delegate				*m_delegate{};
			dispatch_ptr <dispatch_queue_t>		m_worker_queue{};
			dispatch_ptr <dispatch_semaphore_t>	m_process_sema{};
			variant_set							m_factory;
			variant_set							m_prepared_variants;
			variant_vector						m_buffer;
			std::size_t							m_previous_pos{};
			
			data() = default;
			
			data(
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
		};
		
	protected:
		std::mutex							m_buffer_mutex{};
		data								m_d{};
		
	protected:
		void return_node_to_buffer(variant_set::node_type &&node);
		bool get_node_from_buffer(variant_set::node_type &node);

	public:
		variant_buffer() = default;
		
		variant_buffer(variant_buffer &&other)
		{
			using std::swap;
			swap(*this, other);
		}
		
		variant_buffer &operator=(variant_buffer &&other) &
		{
			using std::swap;
			swap(*this, other);
			return *this;
		}
		
		variant_buffer(
			vcf_reader &reader,
			dispatch_ptr <dispatch_queue_t> const &worker_queue,
			variant_buffer_delegate &delegate
		):
			m_d(reader, worker_queue, delegate)
		{
		}
		
		vcf_reader &reader() { return *m_d.m_reader; }
		void read_input();
		void process_input(variant_set &variants);
		void set_delegate(variant_buffer_delegate &delegate) { m_d.m_delegate = &delegate; }
	};
	
	
	inline void swap(variant_buffer &lhs, variant_buffer &rhs) { using std::swap; swap(lhs.m_d, rhs.m_d); }
}

#endif
