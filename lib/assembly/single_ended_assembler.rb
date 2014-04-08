require 'ds'
require 'set'

class Bio::AssemblyGraphAlgorithms::SingleEndedAssembler
  include Bio::FinishM::Logging

  DEFAULT_MAX_TIP_LENGTH = 100

  # Assemble everything in the graph into OrientedNodeTrail objects.
  # Yields an OrientedNodeTrail if a block is
  # given, otherwise returns an array of found paths
  #
  # Options:
  # :min_coverage_of_start_nodes: only start exploring from nodes with this much coverage
  # :min_contig_size: don't print contigs shorter than this (default 500bp)
  # options from #assemble_from are all valid here also
  def assemble(graph, sequences, options={})
    options[:min_contig_size] ||= 500
    options[:max_tip_length] ||= DEFAULT_MAX_TIP_LENGTH

    # Gather a list of nodes to start from
    starting_nodes = []
    if options[:min_coverage_of_start_nodes]
      graph.nodes.each do |node|
        unless node.coverage < options[:min_coverage_of_start_nodes]
          starting_nodes.push node
        end
      end
    else
      starting_nodes = graph.nodes
    end
    log.debug "Found #{starting_nodes.length} nodes to attempt assembly from" if log.debug?

    seen_nodes = Set.new
    # For each starting node, start the assembly process
    dummy_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    starting_nodes.each do |start_node|
      # If we've already covered this node, don't try it again
      next if seen_nodes.include?(start_node.node_id)

      # first find the start of the path
      # Need the start of the path to be not on a short tip,
      # and that it is at the start/end of the path, not the middle.
      # TODO: deal with a fully circular genome where this never ends.
      cross_node = find_nearest_node_on_a_path(node, options[:max_tip_length])
      if cross_node.nil?
        # There's no paths around here, give up.
        seen_nodes << node.node_id
        next
      end
      # Go all the way to the start
      onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
      onode.node = cross_node
      onode.first_side = Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST
      dummy_trail.trail = [onode]
      path, just_visited_onodes = assemble_from(dummy_trail, graph, sequences, options)

      # Then start assembling forwards
      real_start_onode = path[-1]
      real_start_onode.reverse!
      dummy_trail.trail = [real_start_onode]
      path, just_visited_onodes = assemble_from(dummy_trail, graph, sequences, options)

      # Record which nodes have already been visited, so they aren't visited again
      visited_onodes.merge just_visited_onodes.collect{|onode| onode.node_id}

      next if path.length_in_bp < options[:min_contig_size]
      if block_given?
        yield path
      else
        paths.push path
      end
    end
    return paths
  end

  # Given a node, find a closely connected node that is on an assembly path
  # assuming min_contig_size == 0
  def find_nearest_node_on_a_path(node, max_tip_length)
    onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new
    onode.node = node
    onode.first_side = Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST
    raise
  end

  # Assemble considering reads all reads as single ended. Options:
  # :max_tip_length: if a path is shorter than this in bp, then it will be clipped from the path. Default 100
  # :recoherence_kmer: attempt to separate paths by going back to the reads with this larger kmer
  # TODO: deal with a fully circular genome where this never ends.
  def assemble_from(initial_path, graph, sequences, options={})
    options[:max_tip_length] ||= DEFAULT_MAX_TIP_LENGTH

    recoherencer = Bio::AssemblyGraphAlgorithms::SingleCoherentPathsBetweenNodesFinder.new

    path = initial_path.copy
    visited_onodes = []
    dummy_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
    while true
      log.debug "Starting from #{path[-1].to_shorthand}" if log.debug?
      oneighbours = path.neighbours_of_last_node(graph)
      if oneighbours.length == 0
        log.debug "Found a dead end, last node is #{path[-1].to_shorthand}" if log.debug?
        return path, visited_onodes

      elsif oneighbours.length == 1
        to_add = oneighbours[0]
        log.debug "Only one way to go, so going there, to #{to_add.to_shorthand}" if log.debug?
        path.add_oriented_node to_add
        visited_onodes.push to_add

      else
        # Reached a fork (or 3 or 4-fork), which way to go?

        # Remove neighbours that are short tips
        oneighbours.reject! do |oneigh|
          is_tip, visiteds = is_short_tip?(oneigh, graph, options[:max_tip_length])
          if is_tip
            visiteds.each do |onode|
              visited_onodes.push onode
            end
          end
          is_tip
        end

        if oneighbours.length == 0
          log.debug "Found a dead end at a fork, last node is #{path[-1].to_shorthand}" if log.debug?
          return path, visited_onodes
        elsif oneighbours.length == 1
          log.debug "Clipped short tip(s) off, and then there was only one way to go" if log.debug?
          path.add_oriented_node oneighbours[0]
        visited_onodes.push to_add
        elsif options[:recoherence_kmer].nil?
          log.debug "Came across what appears to be a legitimate fork and no recoherence kmer given, so giving up" if log.debug?
          return path, visited_onodes
        else
          unless options[:recoherence_kmer].nil?
            log.debug "Attempting to resolve fork by recoherence" if log.debug?
            oneighbours.select! do |oneigh|
              dummy_trail.trail = path.trail+[oneigh]
              recoherencer.validate_last_node_of_path_by_recoherence(
                dummy_trail,
                options[:recoherence_kmer],
                sequences
                )
            end
          end
          if oneighbours.length == 0
            log.debug "no neighbours passed recoherence, giving up" if log.debug?
            return path
          elsif oneighbours.length == 1
            log.debug "After recoherence there's only one way to go, going there"
            path.add_oriented_node oneighbours[0]
            visited_onodes.push to_add
          else
            log.debug "Still forked after recoherence, so seems to be a legitimate fork, giving up" if log.debug?
            return path, visited_onodes
          end
        end
      end
    end

    return path, visited_onodes
  end

  # Returns false iff there is a path longer than max_tip_length
  # starting at the given oriented_node. Currently works like Dijkstra's
  # shortest path algorithm except that it finds the longest path, not the
  # shortest.
  def is_short_tip?(oriented_node, graph, max_tip_length)
    stack = DS::Stack.new
    first = MaxDistancedOrientedNode.new
    first.onode = oriented_node
    first.distance = oriented_node.node.length_alone
    stack.push first

    cache = {}

    while current_max_distanced_onode = stack.pop
      return false if current_max_distanced_onode.distance > max_tip_length

      current_max_distanced_onode.onode.next_neighbours(graph).each do |oneigh|
        neighbour_distance = current_max_distanced_onode.distance + oneigh.node.length_alone
        next if cache[oneigh.to_settable] and cache[oneigh.to_settable] >= neighbour_distance
        distanced_node = MaxDistancedOrientedNode.new
        distanced_node.onode = oneigh
        distanced_node.distance = neighbour_distance
        cache[oneigh.to_settable] = neighbour_distance
        stack.push distanced_node
      end
    end

    return true, cache.collect{|donode| onode}
  end

  class MaxDistancedOrientedNode
    attr_accessor :onode, :distance
  end
end
