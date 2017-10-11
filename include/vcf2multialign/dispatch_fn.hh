/*
 * Copyright (c) 2016-2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_DISPATCH_FN_HH
#define VCF2MULTIALIGN_DISPATCH_FN_HH

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <dispatch/dispatch.h>
#include <iostream>
#include <string>
#include <stdexcept>


namespace vcf2multialign { namespace detail {
	
	// A virtual destructor is required in dispatch_source_set_event_handler_fn.
	class dispatch_fn_context_base
	{
	public:
		virtual ~dispatch_fn_context_base() {}
	};
	
	template <typename Fn>
	class dispatch_fn_context : public dispatch_fn_context_base
	{
	public:
		typedef Fn function_type;
		
	protected:
		function_type m_fn;
	
	protected:
		static inline void do_call_fn(dispatch_fn_context &ctx)
		{
			try
			{
				ctx.m_fn();
			}
			catch (std::exception const &exc)
			{
				std::cerr << "Caught exception: " << exc.what() << std::endl;
			}
			catch (...)
			{
				std::cerr << "Caught non-std::exception." << std::endl;
			}
		}
	
	public:
		dispatch_fn_context(Fn &&fn):
			m_fn(std::move(fn))
		{
		}
		
		virtual ~dispatch_fn_context() {}
		
		static void cleanup(void *dispatch_context)
		{
			assert(dispatch_context);
			auto *ctx(reinterpret_cast <dispatch_fn_context *>(dispatch_context));
			delete ctx;
		}
		
		static void call_fn_no_delete(void *dispatch_context)
		{
			assert(dispatch_context);
			auto *ctx(reinterpret_cast <dispatch_fn_context *>(dispatch_context));
			do_call_fn(*ctx);
		}
		
		static void call_fn(void *dispatch_context)
		{
			assert(dispatch_context);
			auto *ctx(reinterpret_cast <dispatch_fn_context *>(dispatch_context));
			do_call_fn(*ctx);
			delete ctx;
		}
	};
	
	
	template <typename t_owner, void(t_owner::*t_fn)()>
	void call_member_function(void *ctx)
	{
		t_owner *owner(static_cast <t_owner *>(ctx));
		(owner->*t_fn)();
	}
}}


namespace vcf2multialign {
	// Allow passing pointer-to-member-function to dispatch_async_f without std::function.
	template <typename t_owner>
	class dispatch_caller
	{
	protected:
		t_owner	*m_owner{nullptr};
		
	public:
		dispatch_caller(t_owner *owner): m_owner(owner) { assert(m_owner); }
		
		template <void(t_owner::*t_fn)()>
		void async(dispatch_queue_t queue)
		{
			dispatch_async_f(queue, m_owner, detail::call_member_function <t_owner, t_fn>);
		}
		
		template <void(t_owner::*t_fn)()>
		void barrier_async(dispatch_queue_t queue)
		{
			dispatch_barrier_async_f(queue, m_owner, detail::call_member_function <t_owner, t_fn>);
		}
		
		template <void(t_owner::*t_fn)()>
		void source_set_event_handler(dispatch_source_t source)
		{
			dispatch_set_context(source, m_owner);
			dispatch_source_set_event_handler_f(source, detail::call_member_function <t_owner, t_fn>);
		}
	};
	
	template <typename t_owner>
	dispatch_caller <t_owner> dispatch(t_owner *owner)
	{
		return dispatch_caller <t_owner>(owner);
	}
	
	template <typename Fn>
	void dispatch_async_fn(dispatch_queue_t queue, Fn fn)
	{
		// A new expression doesn't leak memory if the object construction throws an exception.
		typedef detail::dispatch_fn_context <Fn> context_type;
		auto *ctx(new context_type(std::move(fn)));
		dispatch_async_f(queue, ctx, &context_type::call_fn);
	}
	
	template <typename Fn>
	void dispatch_barrier_async_fn(dispatch_queue_t queue, Fn fn)
	{
		typedef detail::dispatch_fn_context <Fn> context_type;
		auto *ctx(new context_type(std::move(fn)));
		dispatch_barrier_async_f(queue, ctx, &context_type::call_fn);
	}
	
	template <typename Fn>
	void dispatch_sync_fn(dispatch_queue_t queue, Fn fn)
	{
		typedef detail::dispatch_fn_context <Fn> context_type;
		auto *ctx(new context_type(std::move(fn)));
		dispatch_sync_f(queue, ctx, &context_type::call_fn);
	}
	
	template <typename Fn>
	void dispatch_source_set_event_handler_fn(dispatch_source_t source, Fn fn)
	{
		typedef detail::dispatch_fn_context <Fn> context_type;
		
		// If the source has been cancelled, dispatch_get_context will return a dangling pointer.
		assert(!dispatch_source_testcancel(source));
		
		{
			// If there is an old context, deallocate it.
			auto *dispatch_context(dispatch_get_context(source));
			if (dispatch_context)
			{
				auto *ctx(reinterpret_cast <detail::dispatch_fn_context_base *>(dispatch_context));
				delete ctx;
			}
		}
		
		auto *new_ctx(new context_type(std::move(fn)));
		dispatch_set_context(source, new_ctx);
		dispatch_source_set_event_handler_f(source, &context_type::call_fn_no_delete);
		dispatch_source_set_cancel_handler_f(source, &context_type::cleanup);
	}
	
	
	template <typename t_dispatch>
	class dispatch_ptr
	{
	protected:
		t_dispatch	m_ptr{};
		
	public:
		dispatch_ptr()
		{
		}
		
		dispatch_ptr(t_dispatch ptr, bool retain = false):
			m_ptr(ptr)
		{
			if (m_ptr && retain)
				dispatch_retain(m_ptr);
		}
		
		~dispatch_ptr()
		{
			if (m_ptr)
				dispatch_release(m_ptr);
		}
		
		dispatch_ptr(dispatch_ptr const &other):
			dispatch_ptr(other.m_ptr, true)
		{
		}
		
		dispatch_ptr(dispatch_ptr &&other):
			m_ptr(other.m_ptr)
		{
			other.m_ptr = nullptr;
		}
		
		bool operator==(dispatch_ptr const &other) const
		{
			return m_ptr == other.m_ptr;
		}
		
		bool operator!=(dispatch_ptr const &other) const
		{
			return !(*this == other);
		}
		
		dispatch_ptr &operator=(dispatch_ptr const &other) &
		{
			if (*this != other)
			{
				if (m_ptr)
					dispatch_release(m_ptr);
				
				m_ptr = other.m_ptr;
				dispatch_retain(m_ptr);
			}
			return *this;
		}
		
		dispatch_ptr &operator=(dispatch_ptr &&other) &
		{
			if (*this != other)
			{
				m_ptr = other.m_ptr;
				other.m_ptr = nullptr;
			}
			return *this;
		}
		
		t_dispatch operator*() { return m_ptr; }
	};
	
	
	template <typename t_dispatch>
	void swap(dispatch_ptr <t_dispatch> &lhs, dispatch_ptr <t_dispatch> &rhs)
	{
		using std::swap;
		swap(lhs.m_ptr, rhs.m_ptr);
	}
}

#endif
