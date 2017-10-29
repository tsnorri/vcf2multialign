/*
 * Copyright (c) 2017 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#ifndef VCF2MULTIALIGN_TASKS_TASK_HH
#define VCF2MULTIALIGN_TASKS_TASK_HH


namespace vcf2multialign {
	
	class task
	{
	public:
		virtual ~task() {}
		virtual void execute() = 0;
	};
}

#endif
