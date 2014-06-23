require 'set'
require 'ds'

# Like DS::PriorityQueue except give the ability to define how priority is given
class DS::AnyPriorityQueue < DS::PriorityQueue
  #Create new priority queue. Internaly uses heap to store elements.
  def initialize
    @store = DS::BinaryHeap.new {|parent,child| yield parent.key, child.key}
  end
end



class Bio::AssemblyGraphAlgorithms::BubblyAssembler < Bio::AssemblyGraphAlgorithms::SingleEndedAssembler
  include Bio::FinishM::Logging

  # Starting at a node within a graph, walk through the graph
  # accepting forks, so long as the fork paths converge within some finite
  # length in the graph (the leash length, measured in number of base pairs).
  #
  # Return an Array of Path arrays, a MetaPath, where each path array are the different paths
  # that can be taken at each fork point
  def assemble_from_node(starting_path, leash_length)
    current_bubble = nil

    metapath = MetaPath.new
    starting_path.each do |oriented_node|
      log.debug "adding onode at the start: #{oriented_node.to_shorthand}" if log.debug?
      metapath << oriented_node
    end

    # Priority queue to determine which path in the bubble to explore next
    # Prioritise things that have lesser numbers (of bp), not greater numbers as is default
    create_new_queue = lambda {DS::AnyPriorityQueue.new {|a,b| a<=b}}
    queue = create_new_queue.call

    while true
      if current_bubble.nil?
        log.debug "Starting a non-bubble from #{metapath.to_shorthand}" if log.debug?
        while oriented_neighbours = metapath.last.next_neighbours(@graph)
          if oriented_neighbours.empty?
            # This is just a straight out dead end, and we can go no further.
            metapath.fate = MetaPath::DEAD_END_FATE
            return metapath
          elsif oriented_neighbours.length == 1
            # Linear thing here, just keep moving forward
            neighbour = oriented_neighbours[0]

            # Stop if a circuit is detected
            if metapath.includes_node?(neighbour)
              metapath.fate = MetaPath::CIRCUIT_FATE
              return metapath
            else
              metapath << neighbour
            end

          else
            # Reached a fork in the graph here, the point of this algorithm, really.
            current_bubble = Bubble.new
            oriented_neighbours.each_with_index do |oneigh, i|
              # Stop if a circuit is detected
              if metapath.includes_node?(oneigh)
                metapath.fate = MetaPath::CIRCUIT_FATE
                return metapath
              end
              current_bubble.add_oriented_node i, oneigh
              queue.enqueue [i, oneigh.node_id, oneigh.first_side], 0
            end
            break #break out of linear path mode while loop
          end
        end


      else
        # We are in a bubble. Go get some.
        path_index, node_id, first_side = queue.shift
        log.debug "In a bubble: #{current_bubble.to_shorthand}" if log.debug?
        log.debug "dequeued path #{path_index}, node #{node_id}/#{first_side}" if log.debug?
        if path_index.nil?
          # finished, we can go no further. Do we treat this as fluff or as the end of a contig?
          raise "hmm"
        else
          # are we now converged?
          if current_bubble.node_converges_bubble?(node_id, first_side)
            # convergement!
            clean_bubble = current_bubble.converge_on(node_id, first_side)
            # convert clean_bubble int an array of OrientedNodeTrail objects
            trails = []
            clean_bubble.each do |trail_of_settabled|
              trail = Bio::Velvet::Graph::OrientedNodeTrail.new
              trail.add_setabled_nodes(trail_of_settabled, @graph)
              trails.push trail
            end
            metapath << trails
            onode = Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(
              @graph.nodes[node_id], first_side
            )
            metapath << onode
            queue = create_new_queue.call
            current_bubble = nil

          else
            # we are not converged. Onwards.
            if !leash_length.nil? and current_bubble.path_length(path_index) > leash_length
              # This bubble did not converge within the specified leash length. So this is the end for this bubble,
              # it must be split up.
              metapath.fate = MetaPath::DIVERGES_FATE
              return metapath
            else
              # Trekking onwards in the bubble. Add this to the path, and queue up the next neighbours
              neighbours = @graph.neighbours_of_node_id(node_id, first_side)
              if neighbours.empty?
                # Dead end. fail? Fluff?
                raise "hmm2"
              else
                # add to the current path, or diverge even further
                new_path_indices = current_bubble.add_one_or_more_neighbours(path_index, neighbours)

                # Add the next neighbours as those to be processed
                new_path_indices.each_with_index do |path_id, i|
                  queue.enqueue(
                    [path_id, neighbours[i].node_id, neighbours[i].first_side],
                    current_bubble.path_length(path_id, @graph)
                    )
                end
              end
            end
          end
        end
      end
    end
  end
end




class Bio::AssemblyGraphAlgorithms::BubblyAssembler::MetaPath < Array
  DIVERGES_FATE = 'diverges'
  DEAD_END_FATE = 'dead end'
  CIRCUIT_FATE = 'circuit'

  # How does this metapath end?
  attr_accessor :fate

  def initialize
    @all_nodes = Set.new
  end

  def last
    self[-1]
  end

  # Returns true if this oriented node is in any path
  def includes_node?(oriented_node)
    @all_nodes.include?(oriented_node.to_settable)
  end

  def <<(oriented_node_or_path_array)
    if oriented_node_or_path_array.kind_of?(Array)
      #its a bubble, an Array of arrays
      oriented_node_or_path_array.each do |trail|
        trail.each do |onode|
          @all_nodes << onode.to_settable
        end
      end
    else
      #its a plain old oriented node
      @all_nodes << oriented_node_or_path_array.to_settable
    end
    super
  end
  alias_method :push, :<<

  def to_shorthand
    to_return = []
    each do |e|
      if e.kind_of?(Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode)
        to_return.push e.to_shorthand
      else
        # e.g.
        to_return.push "(#{e.collect{|e2| e2.to_shorthand}.join(',')})"
      end
    end
    return to_return.join(',')
  end
end

class Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble < Array
  include Bio::FinishM::Logging

  # Return true if the given node is present in all paths, else false
  def node_converges_bubble?(node_id, first_side)
    raise "programming error" if empty?
    to_test = [node_id, first_side]
    each do |other_path|
      unless other_path.key?(to_test)
        log.debug "#{to_test} is not a convergent node because it fails on #{other_path}" if log.debug?
        return false
      end
    end
    log.debug "#{to_test} was a bubble convergent node" if log.debug?
    return true
  end

  # Add an oriented node to the specified path number.
  def add_oriented_node(path_id, oriented_node)
    self[path_id] ||= {}
    # yey ordered hashes in ruby. Assumes nodes are never removed from paths.
    self[path_id][oriented_node.to_settable] = self[path_id].length
    return nil
  end

  # Given a path ID and an Enumerable of oriented nodes, add these to
  # the path in the bubble, splitting it if necessary. Return an Array
  # of path indices corresponding to the path indices that the neighbours
  # were added to.
  #TODO: make a ceiling on the maximal number of nodes allowable in a graph, not just leash length
  def add_one_or_more_neighbours(path_id, neighbours)
    path = self[path_id]
    first = nil
    path_indices = [path_id]
    neighbours.each_with_index do |onode, i|
      if i==0 # no need to copy the first path, can re-use it.
        first = onode
      else
        # Duplicate the path
        new_path = path.to_a.to_h
        new_path_id = length
        self[new_path_id] = new_path
        add_oriented_node(new_path_id, onode)

        # Record which path we are on about here
        path_indices.push i
      end
    end
    # Add the node to the first path
    add_oriented_node(path_id, first)

    return path_indices
  end

  # Return the length of the path in base pairs.
  def path_length(path_id, graph)
    to_return = 0
    self[path_id].collect do |onode_settabled, i|
      to_return += graph.nodes[onode_settabled[0]].length_alone
    end
    to_return
  end

  # The given oriented node is in all paths. Return an array of paths
  # (an array of arrays) that includes all nodes up until the given node
  # (but not including it)
  def converge_on(node_id, first_side)
    to_test = [node_id, first_side]
    collect do |hash|
      new_array = []
      hash.each do |onode_settabled, i|
        if onode_settabled == to_test
          break
        else
          new_array << onode_settabled
        end
      end
      new_array
    end
  end

  def to_shorthand
    to_return = []
    each do |path|
      path_shorthand = []
      path.each do |settabled, i|
        first_side = settabled[1]
        if first_side == :start_is_first
          path_shorthand.push "#{settabled[0]}s"
        else
          path_shorthand.push "#{settabled[0]}e"
        end
      end
      to_return.push path_shorthand.join(',')
    end
    return to_return.join('|')
  end
end
