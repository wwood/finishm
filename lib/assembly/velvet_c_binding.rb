require 'ffi'

class Bio::FinishM::VelvetCBinding
  extend FFI::Library
  include Bio::FinishM::Logging

  ffi_lib File.join(File.dirname(__FILE__),'..','external','libfinishm.so.1.0')

  ############## ProbeNodeFinding ##################################
  # IDnum* extract_best_probe_reads(Graph* graph, IDnum* probeReadIDs, IDnum numProbeReads);
  attach_function :extract_best_probe_reads, [:pointer, :pointer, :int32], :pointer

  ############## ReadToNode ##################################
  #   typedef struct {
  #   IDnum read_id;
  #   IDnum node_id;
  # } ReadIdNodeId; //this is a doubly used structure, once for creation of the data structure, and also in the result of the indexation. This is a little bad, but eh.
  class ReadIdNodeIdStruct < FFI::Struct
    layout :read_id, :int32,
    :node_id, :int32
  end

  # typedef struct {
  #   IDnum num_contents;
  #   IDnum num_reads;
  #   IDnum* index;
  #   ReadIdNodeId* contents;
  # } ReadIdToNodeIdLookupTable;
  class ReadIdToNodeIdLookupTableStruct < FFI::Struct
    layout :num_contents, :int32,
    :num_reads, :int32,
    :index, :pointer,
    :contents, :pointer
  end

  # typedef struct {
  #   IDnum num_nodes;
  #   ReadIdNodeId* read_ids_node_ids;
  # } ReadIdToNodeIdIndexation;
  class ReadIdToNodeIdIndexationStruct < FFI::Struct
    layout :num_nodes, :int32,
    :read_ids_node_ids, ReadIdNodeIdStruct.ptr
  end

  # ReadIdToNodeIdIndexation getReadIdToNodeIdIndexation(ReadIdToNodeIdLookupTable* lookupTable, IDnum readID);
  attach_function :getReadIdToNodeIdIndexation, [:pointer, :int32], ReadIdToNodeIdIndexationStruct.by_value
  alias_method :get_read_id_to_node_id_indexation, :getReadIdToNodeIdIndexation

  # ReadIdToNodeIdLookupTable* readReadIdToNodeIdLookupTable(char* fileName);
  attach_function :readReadIdToNodeIdLookupTable, [:string], :pointer
  alias_method :read_read_id_to_node_id_lookup_table, :readReadIdToNodeIdLookupTable
end


