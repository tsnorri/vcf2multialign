/*
 * Copyright (c) 2021 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_UTILITY_DISPATCH_EXIT_GUARD_HH
#define VCF2MULTIALIGN_UTILITY_DISPATCH_EXIT_GUARD_HH

#include <libbio/dispatch.hh>


namespace vcf2multialign {
	
	// May be inherited or used as a data member for sending std::exit to the main queue.
	struct dispatch_exit_guard
	{
		dispatch_exit_guard() = default;
		// Prevent copying and moving.
		dispatch_exit_guard(dispatch_exit_guard const &) = delete;
		dispatch_exit_guard(dispatch_exit_guard &&) = delete;
		dispatch_exit_guard &operator=(dispatch_exit_guard const &) = delete;
		dispatch_exit_guard &operator=(dispatch_exit_guard &&) = delete;
		
		~dispatch_exit_guard()
		{
			dispatch_async(dispatch_get_main_queue(), ^{
				std::exit(EXIT_SUCCESS);
			});
		}
	};
	
	
	// Contains the correct data member order s.t. dispatch_exit_guard will be destroyed last.
	template <typename t_type>
	struct dispatch_exit_guard_helper
	{
		dispatch_exit_guard	guard; // Will be deallocated last.
		t_type				value;
		
		dispatch_exit_guard_helper() = default;
		
		template <typename ... t_args>
		dispatch_exit_guard_helper(t_args && ... args):
			guard(),
			value(std::forward <t_args>(args)...)
		{
		}
	};
}

#endif
