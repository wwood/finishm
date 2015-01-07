class Bio::FinishM::ReadToNode
  include Bio::FinishM::Logging

  def initialize(filename)
    @bindings = Bio::FinishM::VelvetCBinding.new
    log.debug "Reading ReadToNode file #{filename}.."
    raise "Unable to find readToNode binary file" unless File.exist?(filename)
    @read_to_node = @bindings.read_read_id_to_node_id_lookup_table(filename)
    log.debug "Finished reading ReadToNode file"
    @cache = {}
  end

  # Return an array of node IDs that include the given read id
  def [](read_id)
    # cache
    cache = @cache[read_id]
    return cache unless cache.nil?

    res = @bindings.get_read_id_to_node_id_indexation(@read_to_node, read_id)

    #   # typedef struct {
    #   #   IDnum num_nodes;
    #   #   ReadIdNodeId* read_ids_node_ids;
    #   # } ReadIdToNodeIdIndexation;
    #   class ReadIdToNodeIdIndexationStruct < FFI::Struct
    #     layout :num_nodes, :int32,
    #     :read_ids_node_ids, :pointer
    #   end
    to_return = []
    structs = FFI::Pointer.new(Bio::FinishM::VelvetCBinding::ReadIdNodeIdStruct, res[:read_ids_node_ids].pointer)
    0.upto(res[:num_nodes]-1) do |i|
      to_return << Bio::FinishM::VelvetCBinding::ReadIdNodeIdStruct.new(structs[i])[:node_id].abs
    end
    @cache[read_id] = to_return
    return to_return
  end
end
