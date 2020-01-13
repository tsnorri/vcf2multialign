/*
 * Copyright (c) 2019 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <stdexcept>
#include "msa_combiner.hh"


%% machine alignment_parser;
%% write data;


namespace vcf2multialign {
	
	void msa_combiner::prepare_msa_parser()
	{
		%% write init;
	}
	
	
	void msa_combiner::parse_msa(aligned_character_pack const &pack)
	{
		%%{
			variable cs		m_fsm.cs;
			variable p		m_fsm.p;
			variable pe		m_fsm.pe;
			variable eof	m_fsm.pe;
			
			nt			= [ACGTN];	# Add [acgtn] to allow lowecase.
			ntg			= nt | '-';
			
			a_			= [A];		# Add [a] etc. to allow lowercase.
			c_			= [C];
			g_			= [G];
			t_			= [T];
			n_			= [N];
			not_a		= nt - a_;
			not_c		= nt - c_;
			not_g		= nt - g_;
			not_t		= nt - t_;
			not_n		= nt - n_;
			
			same		= /AA/ | /CC/ | /GG/ | /TT/ | /NN/;	# Add /i to allow lowercase.
			diff		= (a_ . not_a) | (c_ . not_c) | (g_ . not_g) | (t_ . not_t) | (n_ . not_n);
			both_nt		= nt{2};				# Any two non-gap.
			
			action deletion_continue_r {
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_deletion;
				fbreak;
			}
			
			action main_new_deletion_r {
				m_current_segment.reset(pack, segment_type::DELETION);
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_deletion;
				fbreak;
			}
			
			action main_new_mismatch_b {
				m_current_segment.reset(pack, segment_type::MISMATCH);
				m_current_segment.ref.string += pack.ref.character;
				m_current_segment.alt.string += pack.alt.character;
				fnext handle_mismatch;
				fbreak;
			}
			
			action main_new_mixed_b {
				m_current_segment.reset(pack, segment_type::MIXED);
				m_current_segment.ref.string += pack.ref.character;
				m_current_segment.alt.string += pack.alt.character;
				fnext handle_mixed;
				fbreak;
			}
			
			action main_new_mixed_r {
				m_current_segment.reset(pack, segment_type::MIXED);
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_mixed;
				fbreak;
			}
			
			action main_new_match_b {
				m_current_segment.reset(pack, segment_type::MATCH);
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_match;
				fbreak;
			}
			
			action match_continue_b {
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_match;
				fbreak;
			}
			
			action match_continue_g {
				fnext handle_match;
				fbreak;
			}
			
			action mismatch_continue_b {
				m_current_segment.ref.string += pack.ref.character;
				m_current_segment.alt.string += pack.alt.character;
				fnext handle_mismatch;
				fbreak;
			}
			
			action mismatch_continue_g {
				fnext handle_mismatch;
				fbreak;
			}
			
			action mismatch_mark_mixed_a {
				m_current_segment.type = segment_type::MIXED;
				m_current_segment.alt.string += pack.alt.character;
				fnext handle_mixed;
				fbreak;
			}
			
			action mixed_continue_a {
				m_current_segment.alt.string += pack.alt.character;
				fnext handle_mixed;
				fbreak;
			}
			
			action mixed_continue_g {
				fnext handle_mixed;
				fbreak;
			}
			
			action mixed_continue_r {
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_mixed;
				fbreak;
			}
			
			action start_new_del_r {
				push_current_segment();
				m_current_segment.reset(pack, segment_type::DELETION);
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_deletion;
				fbreak;
			}
			
			action start_new_match_b {
				push_current_segment();
				m_current_segment.reset(pack, segment_type::MATCH);
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_match;
				fbreak;
			}
			
			action start_new_mismatch_b {
				push_current_segment();
				m_current_segment.reset(pack, segment_type::MISMATCH);
				m_current_segment.ref.string += pack.ref.character;
				m_current_segment.alt.string += pack.alt.character;
				fnext handle_mismatch;
				fbreak;
			}
			
			action start_new_mixed_b {
				push_current_segment();
				m_current_segment.reset(pack, segment_type::MIXED);
				m_current_segment.ref.string += pack.ref.character;
				m_current_segment.alt.string += pack.alt.character;
				fnext handle_mixed;
				fbreak;
			}
			
			action start_new_mixed_r {
				push_current_segment();
				m_current_segment.reset(pack, segment_type::MIXED);
				m_current_segment.ref.string += pack.ref.character;
				fnext handle_mixed;
				fbreak;
			}
			
			action handle_error {
				throw std::runtime_error("Unexpected character");
			}
			
			action handle_unexpected_state {
				throw std::runtime_error("Unexpected state");
			}
			
			# The following machine definitions should handle all combinations of four characters.
			# “Same” or “matching” segment refers to the property that for every alt non-gap character there is a matching ref character.
			# Otherwise the segment is considered non-matching. In some cases, this property cannot be determined by looking ahead just one character.
			
			# Ref	X-	X-	X-	X-
			# Alt	X-	X+	Y*	-*
			same_next_both_g		= (same . '--');
			same_next_ref_g			= (same . '-' . nt);
			diff_next_ref_g			= (diff . '-' . ntg);
			alt_g_next_ref_g		= (nt . '--' . ntg);
			
			# Ref	X+
			# Alt	-*
			alt_g_next_ref_c		= (nt . '-' . nt . ntg);
			
			# Ref	X+	X+
			# Alt	X*	Y*
			# Non-gap alt characters consistently match / do not match ref.
			same_next_ref_nt		= (same . nt . ntg);
			diff_next_ref_nt		= (diff . nt . ntg);
			
			# Ref	-	-
			# Alt	-	X
			# ref_g may come after same_next_ref_g, diff_next_ref_g, alt_g_next_ref_g.
			both_g					= ('--');
			ref_g					= ('-' . nt);
			
			# Make a static assertion fail if there are any character combinations that were not listed above.
			action report_unhandled {
				static_assert(false, "Found a character combination that was not handled.");
			}
			unhandled = ((ntg{4}) - (
				same_next_both_g		|
				same_next_ref_g			|
				diff_next_ref_g			|
				alt_g_next_ref_g		|
				
				alt_g_next_ref_c		|
				
				same_next_ref_nt		|
				diff_next_ref_nt		|
				
				(both_g . ntg{2})		|	# Corresponds to both_g in the machines below but would match more without the ntg{2}.
				(ref_g . ntg{2})			# Corresponds to ref_g in the machines below but would match more without the ntg{2}.
			));
			fail_on_unhandled_combinations := unhandled $(report_unhandled);
		
			# In the following cases start a new segment at the marked position b.c.
			# the following characters cannot be known and this reduces the number
			# of special cases a bit (1).
			# Ref	ACG-	ACG-
			# Alt	ACG-	TTG-
			#		  ^		  ^
			# Action suffixes: _b: both have non-gap, _r: ref has non-gap, _a: alt has non-gap, _g: both have gap.
			
			common = (
				# Since the next ref character is a gap, a new segment needs to be created since we don’t know the subsequent characters.
				(same_next_both_g	@(start_new_mixed_b))		|	# We currently don’t change a matching segment to mixed, so mark the segment for post-processing here.
				(same_next_ref_g	@(start_new_mixed_b))		|
				(diff_next_ref_g	@(start_new_mismatch_b))	|
				(alt_g_next_ref_g	@(start_new_mixed_r))
			);
			
			handle_match := common | (
				(alt_g_next_ref_c	@(start_new_del_r))			|	# Deletion in alt w.r.t. ref.
				
				(same_next_ref_nt	@(match_continue_b))		|
				(diff_next_ref_nt	@(start_new_mismatch_b))	|
				
				(both_g				@(handle_unexpected_state))	|
				(ref_g				@(handle_unexpected_state))
			) $err(handle_error);
			
			handle_mismatch := common | (
				(alt_g_next_ref_c	@(start_new_del_r))			|	# Deletion in alt w.r.t. ref.
	
				(same_next_ref_nt	@(start_new_match_b))		|
				(diff_next_ref_nt	@(mismatch_continue_b))		|
	
				(both_g				@(mismatch_continue_g))		|
				(ref_g				@(mismatch_mark_mixed_a))
			) $err(handle_error);
			
			handle_mixed := common | (
				(alt_g_next_ref_c	@(start_new_del_r))			|	# Deletion in alt w.r.t. ref.

				(same_next_ref_nt	@(start_new_match_b))		|
				(diff_next_ref_nt	@(start_new_mismatch_b))	|
				
				(both_g				@(mixed_continue_g))		|
				(ref_g				@(mixed_continue_a))
			) $err(handle_error);
			
			handle_deletion := common | (
				(alt_g_next_ref_c	@(deletion_continue_r))		|
				
				(same_next_ref_nt	@(start_new_match_b))		|
				(diff_next_ref_nt	@(start_new_mismatch_b))	|
				
				(both_g				@(handle_unexpected_state))	|
				(ref_g				@(handle_unexpected_state))
			) $err(handle_error);
			
			# Starting machines. Ref is guaranteed to be non-gap. See common.
			# Ref	X-	X-	X-	X-	X+	X+	X+
			# Alt	X-	X+	Y*	-*	-*	X*	Y*
			#same_next_both_g
			#same_next_ref_g
			#diff_next_ref_g
			#alt_g_next_ref_g
			#alt_g_next_ref_c
			#same_next_ref_nt
			#diff_next_ref_nt
			main_unexpected = ntg{4} - (same_next_both_g | same_next_ref_g | diff_next_ref_g | alt_g_next_ref_g | alt_g_next_ref_c | same_next_ref_nt | diff_next_ref_nt);
			
			main := (
				(same_next_both_g	@(main_new_match_b))		|
				(same_next_ref_g	@(main_new_mixed_b))		|
				(diff_next_ref_g	@(main_new_mismatch_b))		|
				(alt_g_next_ref_g	@(main_new_mixed_r))		|
				(alt_g_next_ref_c	@(main_new_deletion_r))		|
				(same_next_ref_nt	@(main_new_match_b))		|
				(diff_next_ref_nt	@(main_new_mismatch_b))		|
				(main_unexpected	@(handle_unexpected_state))
			) $err(handle_error);
			
			write exec noend;
		}%%
	}
}
