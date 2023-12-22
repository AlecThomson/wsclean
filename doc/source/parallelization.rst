Parallization options
=====================

WSClean can exploit parallelism by using multiple threads and by using multiple
MPI processes. The MPI processes may run at different compute nodes, which
enables distributed processing.

Multi-threading
---------------
WSClean has several command-line options for configuring multi-threading:

* Number of cores to use: ``-j <N>``

  Tells WSClean to use ``N`` threads. By default, WSClean uses as many threads
  as there are cores in the system.
  WSClean uses this setting everywhere possible, including
  when WSClean uses external libraries, such as IDG and Radler.

* Do :doc:`parallel deconvolution<parallel_deconvolution>`:
  ``-parallel-deconvolution <maxsize>``

  Tells WSClean to do parallel deconvolution where ``maxsize`` is the maximum
  subimage size. For more information, see
  :doc:`parallel deconvolution<parallel_deconvolution>`.

* Reduce memory usage for deconvolution: ``-deconvolution-threads <N>``

  On machines with a large number of cores, using less deconvolution threads
  helps if WSClean runs out of memory during deconvolution. WSClean allocates
  memory for each thread, e.g., for storing the sub-image a number of times.

  This option tells WSClean to use a maximum of ``N`` threads during
  deconvolution. The default value is the number of threads as specified using
  the ``-j`` argument.

* Enable parallel reordering: ``-parallel-reordering <N>``

  Tells WSClean to use parallelism when reordering
  the input measurement sets on disk.
  WSClean uses one reordering task for each input measurements set, and
  executes up to ``N`` tasks in parallel.

  Reordering is bound by disk speed when the number of cores is high.
  Using 4 reordering threads is generally a good balance between disk and CPU
  speed. By default, WSClean therefore uses 4 reordering threads, even if ``-j``
  is specified. The optimal number of reordering threads depends on disk speed.

* Enable parallel gridding: ``-parallel-gridding <N>``

  Tells WSClean to use ``N`` threads during gridding and degridding.
  By default, parallel gridding is disabled, even if ``-j`` is specified.
  Parallel gridding does not work with the ``idg`` gridder, or when using MPI
  (see below).

  The gridders themselves also exploit parallelism internally, however,
  gridders can't always scale well to all cores, in particular when using beam
  or h5parm corrections. In these cases, using parallel gridding might yield
  better performance.

  All gridders get an amount of threads equal to the total number of threads,
  divided by the number of parallel gridders, rounded up. For example,
  ``-parallel-gridding 3 -j 16`` yields 3 parallel gridders that use
  6 threads each.

Multi-node processing (MPI)
---------------------------

WSClean can be compiled with MPI support. If enabled, the compiler will produce
a ``wsclean-mp`` binary which can exploit parallelism across multiple compute
nodes. This binary supports the same command-line options as
the regular ``wsclean`` binary. The main difference is that it performs
parallel gridding using multiple MPI processes. Parallel gridding using multiple
threads, using the ``-parallel-gridding`` option, is not possible with
``wsclean-mp``.

The other multi-threading arguments (described above) apply to each MPI process.
For example, when using ``-j 4``, each MPI process will use 4 threads, even
if these processes run on the same compute node.

Note that WSClean only uses MPI during gridding. Other parts, such as
deconvolution and reordering, only use the main MPI process and therefore only
use multi-threading.

When using a single compute node, using MPI is discouraged, since multi-threaded
gridding using the ``-parallel-gridding`` is more efficient. When using multiple
nodes, the best performance is normally achieved using one MPI process per node
and a number of threads equal to the number of cores per node.

Processing related facets on one node (Work in progress!)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When faceting is enabled, WSClean normally schedules gridding tasks for
different facets at different compute nodes.
It then uses a global writer lock for synchronizing updating the visibilities.

For improving scalability, WSClean will soon support scheduling the tasks
for all facets for a single output channel at a fixed compute node.
The ``-parallel-gridding`` argument then specifies how many facets each node
should process in parallel, using multiple threads. Instead of the global
writer lock, WSClean then uses a local writer lock on each node.

The ``-channel-to-node`` argument specifies the mapping of output channels to
node. For example, ``-channel-to-node 0,0,1,2`` schedules the tasks for
output channel indices 0, 1, 2, 3 at compute nodes 0, 0, 1, 2, respectively.
The length of the list must equal the number of output channels
(``-channels-out`` argument).
If ``-no-work-on-master`` is specified, the list may not contain ``0``.

By default, WSClean uses a round-robin distribution of output channels to nodes,
because later channels are more expensive to grid.
For example, with 10 output channels, 4 nodes, and an active main node, the
mapping is 0,1,2,3,0,1,2,3,0,1.