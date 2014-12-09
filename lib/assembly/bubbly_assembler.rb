require 'set'
require 'ds'

# Like DS::PriorityQueue except give the ability to define how priority is given
class DS::AnyPriorityQueue < DS::PriorityQueue
  #Create new priority queue. Internaly uses heap to store elements.
  def initialize
    @store = DS::BinaryHeap.new {|parent,child| yield parent.key, child.key}
  end

  def each
    @store.to_a.each do |pair|
      yield pair.value
    end
  end
end



class Bio::AssemblyGraphAlgorithms::BubblyAssembler < Bio::AssemblyGraphAlgorithms::SingleEndedAssembler
  include Bio::FinishM::Logging

  DEFAULT_MAX_BUBBLE_LENGTH = 500
  DEFAULT_BUBBLE_NODE_COUNT_LIMIT = 20 #so, so very 'un-educated' guess

  def initialize(graph, assembly_options={})
    opts = assembly_options
    opts[:max_bubble_length] ||= DEFAULT_MAX_BUBBLE_LENGTH
    opts[:bubble_node_count_limit] ||= DEFAULT_BUBBLE_NODE_COUNT_LIMIT
    super graph, opts
  end

  # Starting at a node within a graph, walk through the graph
  # accepting forks, so long as the fork paths converge within some finite
  # length in the graph (the leash length, measured in number of base pairs).
  #
  # Return an Array of Path arrays, a MetaPath, where each path array are the different paths
  # that can be taken at each fork point
  def assemble_from(starting_path, visited_nodes=Set.new)
    leash_length = @assembly_options[:max_bubble_length]
    current_bubble = nil
    log.info "Assembling from: #{starting_path.to_shorthand}" if log.info?

    metapath = MetaPath.new
    starting_path.each do |oriented_node|
      log.debug "adding onode at the start: #{oriented_node.to_shorthand}" if log.debug?
      metapath << oriented_node
    end

    # Keep track of nodes visited in this trajectory already so circuits can be avoided
    visited_oriented_node_settables = Set.new
    starting_path.each do |e|
      if e.kind_of?(Bubble)
        e.oriented_nodes do |onode|
          visited_oriented_node_settables << onode.to_settable
        end
      else
        visited_oriented_node_settables << e.to_settable
      end
    end

    current_mode = :linear # :linear, :bubble, or :finished

    while current_mode != :finished
      if current_mode == :linear
        log.debug "Starting a non-bubble from #{metapath.to_shorthand}" if log.debug?
        while true
          oriented_neighbours = metapath.last_oriented_node.next_neighbours(@graph)
          log.debug "Found oriented neighbours #{oriented_neighbours.collect{|onode| onode.to_shorthand} }" if log.debug?

          legit_neighbours = nil
          # Cut off tips unless it is the only way
          if oriented_neighbours.length == 1
            legit_neighbours = oriented_neighbours
          else
            legit_neighbours = oriented_neighbours.reject do |oneigh|
              is_tip, visiteds = is_short_tip?(oneigh)
              if is_tip
                visited_oriented_node_settables << oneigh.to_settable
                visiteds.each do |v|
                  visited_oriented_node_settables << v
                end
              end
              log.debug "neighbour #{oneigh.to_shorthand} is_tip? #{is_tip}" if log.debug?
              is_tip
            end
          end

          if legit_neighbours.empty?
            # This is just a straight out dead end, and we can go no further.
            log.debug "Dead end reached" if log.debug?
            metapath.fate = MetaPath::DEAD_END_FATE
            current_mode = :finished
            break
          elsif legit_neighbours.length == 1
            # Linear thing here, just keep moving forward
            neighbour = legit_neighbours[0]

            # Stop if a circuit is detected
            if visited_oriented_node_settables.include?(neighbour.to_settable)
              log.debug "Detected regular circuit by running into #{neighbour.to_settable}" if log.debug?
              metapath.fate = MetaPath::CIRCUIT_FATE
              current_mode = :finished
              break
            else
              visited_oriented_node_settables << neighbour.to_settable
              metapath << neighbour
            end

          else
            # Reached a fork in the graph here, the point of this algorithm, really.
            current_bubble = Bubble.new
            log.debug "Starting a bubble forking from metapath #{metapath.to_shorthand}" if log.debug?

            legit_neighbours.each do |oneigh|
              # Stop if a circuit is detected
              if visited_oriented_node_settables.include?(oneigh.to_settable)
                metapath.fate = MetaPath::CIRCUIT_FATE
                current_mode = :finished
                break
              else
                new_problem = DynamicProgrammingProblem.new
                new_problem.distance = 0
                new_path = Bio::Velvet::Graph::OrientedNodeTrail.new
                new_path.add_oriented_node oneigh
                new_problem.path = new_path
                new_problem.ubiquitous_oriented_nodes = Set.new
                new_problem.ubiquitous_oriented_nodes << oneigh.to_settable

                log.debug "Adding problem to bubble: #{new_problem}" if log.debug?

                current_bubble.enqueue new_problem
                current_mode = :bubble
              end
            end
            break
          end
        end


      elsif current_mode == :bubble
        # We are in a bubble. Go get some.
        log.debug "entering bubble mode" if log.debug?

        # next problem = queue.shift. while distance of next problem is not beyond the leash length
        while current_mode == :bubble
          problem = current_bubble.shift
          log.debug "Dequeued #{problem.to_shorthand}" if log.debug?

          if problem.nil?
            # Getting here seems improbable if not impossible.
            # The current bubble doesn't converge and just has short tips at the end, don't add it to the metapath
            metapath.fate = MetaPath::DEAD_END_FATE
            current_mode = :finished
            log.debug "Reached a dead end, ignoring this path" if log.debug?
            break
          elsif !leash_length.nil? and problem.distance > leash_length
            # The current bubble doesn't converge, don't add it to the metapath
            metapath.fate = MetaPath::DIVERGES_FATE
            current_mode = :finished
            log.debug "Bubble is past the leash length of #{leash_length}, giving up" if log.debug?
            break
          elsif current_bubble.convergent_on?(problem)
            log.debug "Bubble #{current_bubble.to_shorthand} convergent on #{problem.to_shorthand}" if log.debug?
            current_bubble.converge_on problem
            # convergement!
            # Bubble ended in a convergent fashion

            # detect circuits
            if current_bubble.circuitous?
              metapath.fate = MetaPath::CIRCUIT_WITHIN_BUBBLE_FATE
              current_mode = :finished
              break
            else
              metapath << current_bubble
              # Add the nodes in the bubble to the list of visited nodes
              current_bubble.oriented_nodes do |onode|
                visited_oriented_node_settables << onode.to_settable
              end

              current_bubble = nil
              current_mode = :linear
              break
            end
          else
            # otherwise we must search on in the bubble
            # get all neighbours that are not short tips
            log.debug "Bubble not convergent on #{problem.to_shorthand}" if log.debug?
            legit_neighbours = problem.path.neighbours_of_last_node(@graph).reject do |oneigh|
              is_tip, visiteds = is_short_tip?(oneigh)
              if is_tip
                visited_oriented_node_settables << oneigh.to_settable
                visiteds.each do |v|
                  visited_oriented_node_settables << v
                end
              end
              log.debug "neighbour #{oneigh.to_shorthand} is_tip? #{is_tip}" if log.debug?
              is_tip
            end

            if legit_neighbours.length == 0
              # this is a kind of 'long' tip, possibly unlikely to happen much.
              # Forget about it and progress to the next problem having effectively
              # removed it from the bubble
              log.debug "Found no neighbours to re-enqueue" if log.debug?
            else
              legit_neighbours.each do |oneigh|
                new_problem = DynamicProgrammingProblem.new
                new_problem.distance = problem.distance + problem.path[-1].node.length_alone
                new_path = problem.path.copy
                new_path.add_oriented_node oneigh
                new_problem.path = new_path
                new_problem.ubiquitous_oriented_nodes = Set.new problem.ubiquitous_oriented_nodes
                if new_problem.ubiquitous_oriented_nodes.include?(oneigh.to_settable)
                  log.debug "Found circuit within bubble on #{oneigh}, giving up" if log.debug?
                  metapath.fate = MetaPath::CIRCUIT_WITHIN_BUBBLE_FATE
                  current_mode = :finished
                  break
                else
                  new_problem.ubiquitous_oriented_nodes << oneigh.to_settable

                  current_bubble.enqueue new_problem
                  log.debug "Enqueued #{new_problem.to_shorthand}, total nodes now #{current_bubble.num_known_problems}" if log.debug?

                  # check to make sure we aren't going overboard in the bubbly-ness
                  if !@assembly_options[:bubble_node_count_limit].nil? and current_bubble.num_known_problems > @assembly_options[:bubble_node_count_limit]
                    log.debug "Too complex a bubble detected, giving up" if log.debug?
                    metapath.fate = MetaPath::NODE_COUNT_LIMIT_REACHED
                    current_mode = :finished
                    break
                  end
                end
              end
            end
          end
        end
      else
        raise "Programming error: Unexpected mode: #{current_mode}"
      end

      log.debug "Reached end of main loop in mode #{current_mode}" if log.debug?
    end

    return metapath, visited_oriented_node_settables
  end

  def seen_last_in_path?(path, seen_nodes)
    last = path[-1]
    if last.kind_of?(Bubble)
      return remove_seen_nodes_from_end_of_path(path, seen_nodes).length < path.length
    else
      return seen_nodes.include?(path[-1].to_settable)
    end
  end


  def remove_seen_nodes_from_end_of_path(path, seen_nodes)
    log.debug "Removing from the end of the path #{path.to_shorthand} any nodes in set of length #{seen_nodes.length}" if log.debug?

    node_seen = lambda do |oriented_node|
      seen_nodes.include?([oriented_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST]) or
      seen_nodes.include?([oriented_node.node_id, Bio::Velvet::Graph::OrientedNodeTrail::END_IS_FIRST])
    end

    while !path.empty?
      last_node_or_bubble_index = path.length-1
      last_node_or_bubble = path[last_node_or_bubble_index]

      delete = false
      if last_node_or_bubble.kind_of?(Bubble)
        last_node_or_bubble.oriented_nodes do |onode|
          if node_seen.call(onode)
            delete = true
            break
          end
        end
      else
        delete = node_seen.call(last_node_or_bubble)
      end

      if delete
        path.delete_at last_node_or_bubble_index
      else
        # Last node is not previously seen, chop no further.
        break
      end
    end

    return path
  end

  def is_short_tip?(element)
    if element.kind_of?(Bubble)
      return false
    else
      super
    end
  end



  class MetaPath
    DIVERGES_FATE = 'diverges'
    DEAD_END_FATE = 'dead end'
    CIRCUIT_FATE = 'circuit'
    NODE_COUNT_LIMIT_REACHED = 'too many nodes in bubble'
    CIRCUIT_WITHIN_BUBBLE_FATE = 'circuit within bubble'

    # How does this metapath end?
    attr_accessor :fate

    include Enumerable

    def initialize
      @internal_array = []
    end

    def each
      @internal_array.each do |e|
        yield e
      end
    end

    def [](index)
      @internal_array[index]
    end

    def delete_at(index)
      @internal_array.delete_at index
    end

    def empty?
      @internal_array.empty?
    end

    def last_oriented_node
      e = @internal_array[-1]
      if e.kind_of?(Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode)
        return e
      else
        # it is a bubble
        return e.converging_oriented_node
      end
    end

    def <<(oriented_node_or_bubble)
      @internal_array << oriented_node_or_bubble
    end
    alias_method :push, :<<

    def to_shorthand
      @internal_array.collect{|e| e.to_shorthand}.join(',')
    end

    def reverse!
      # Do regular reversal
      @internal_array.reverse!

      # Reverse all the internal parts
      @internal_array.each do |e|
        e.reverse!
      end

      return nil
    end

    def length
      @internal_array.length
    end

    def [](index)
      @internal_array[index]
    end

    def delete_at(index)
      @internal_array.delete_at(index)
    end

    # Yield all oriented nodes anywhere in the regular or bubble
    # bits.
    def each_oriented_node
      @internal_array.each do |e|
        if e.kind_of?(Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble)
          e.oriented_nodes.each do |onode|
            yield onode
          end
        else
          yield e
        end
      end
    end

    def length_in_bp
      sum = 0
      each do |e|
        if e.kind_of?(Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble)
          sum += e.length_in_bp
        else
          sum += e.node.length_alone
        end
      end
      return sum
    end

    def reference_trail
      trail = Bio::Velvet::Graph::OrientedNodeTrail.new

      trail.trail = collect do |e|
        if e.kind_of?(Bio::AssemblyGraphAlgorithms::BubblyAssembler::Bubble)
          e.reference_trail.trail
        else
          e
        end
      end.flatten

      return trail
    end

    def sequence
      reference_trail.sequence
    end

    def coverage
      reference_trail.coverage
    end
  end



  class Bubble
    include Bio::FinishM::Logging

    # The DynamicProgrammingProblem this bubble converges on
    attr_reader :converging_oriented_node_settable, :is_reverse

    def initialize
      @queue = DS::AnyPriorityQueue.new {|a,b| a<=b}
      @known_problems = {}
      @current_problems = Set.new
    end

    # Return the next closest dynamic programming problem,
    # removing it from the bubble
    def shift
      prob = @queue.shift
      @current_problems.delete prob.to_settable
      return prob
    end

    def enqueue(dynamic_programming_problem)
      settable = dynamic_programming_problem.to_settable

      @known_problems[settable] ||= []
      @known_problems[settable].push dynamic_programming_problem

      # Don't re-enqueue the problem if it is already in the queue
      unless @current_problems.include?(settable)
        @queue.enqueue dynamic_programming_problem, dynamic_programming_problem.distance
        @current_problems << settable
      end
    end

    # return true if the given problem converges the bubble, else false
    def convergent_on?(dynamic_programming_problem)
      settable =  dynamic_programming_problem.path[-1].to_settable
      @queue.each do |problem| #convergent until not
        return false unless problem.ubiquitous_oriented_nodes.include?(settable)
      end
      return true
    end

    # Finish off the bubble, assuming convergent_on? the given problem == true
    def converge_on(dynamic_programming_problem)
      @converging_oriented_node_settable = dynamic_programming_problem.to_settable
      #free some memory
      @queue = nil
      @current_problems = nil
    end

    # yield or failing that return an Array of the list of oriented_nodes found
    # in at least one path in this (presumed converged) bubble
    def oriented_nodes
      raise unless converged?
      seen_nodes = {}
      stack = DS::Stack.new
      @known_problems[@converging_oriented_node_settable].each do |initial_solution|
        stack.push initial_solution.path[-1]
      end

      while onode = stack.pop
        setable = onode.to_settable
        next if seen_nodes.key?(setable)

        if block_given?
          if @is_reverse
            yield onode.reverse
          else
            yield onode
          end
        end

        seen_nodes[setable] = onode

        # queue neighbours
        @known_problems[setable].each do |dpp|
          stack.push dpp.path[-2] unless dpp.path.length < 2
        end
      end

      return nil if block_given?
      return seen_nodes.values
    end

    def num_known_problems
      @known_problems.length
    end

    # Iterate over the paths returning each as an OrientedNodeTrail.
    # Assumes the path is convergent.
    def each_path
      raise unless converged?
      problems_to_yield = DS::Queue.new
      @known_problems[@converging_oriented_node_settable].each do |initial_solution|
        problems_to_yield.push [initial_solution.path, []]
        #log.debug "Pushed to stack #{initial_solution.path.to_shorthand}" if log.debug?
      end

      while path_halves = problems_to_yield.dequeue
        direct_node_trail = path_halves[0]
        second_half = path_halves[1]
        #log.debug "Dequeued #{direct_node_trail.to_shorthand} and [#{second_half.collect{|o| o.to_shorthand}.join(',') }]" if log.debug?

        # yield the direct path
        direct_path = Bio::Velvet::Graph::OrientedNodeTrail.new
        direct_path.trail = [direct_node_trail.trail, second_half].flatten
        if @is_reverse
          yield direct_path.reverse
        else
          yield direct_path
        end

        # go down the path, looking for other paths
        unless direct_node_trail.trail.length < 3
          new_head_onode = direct_node_trail.trail[-2]

          new_problems = @known_problems[new_head_onode.to_settable]
          #log.debug "Found new problems: #{new_problems.collect{|prob| prob.to_shorthand}.join(' ') }" if log.debug?
          if new_problems.length > 1
            raise CircuitousPathDetected if second_half.include?(new_head_onode)
            new_second_half = [new_head_onode]+[direct_node_trail.trail[-1]]+second_half
            new_problems.each do |new_problem|
              # Don't enqueue the path that was just printed
              next if direct_path.trail[0...-1] == new_problem.path.trail

              # TODO: deal with circuits
              new_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
              new_trail.trail = new_problem.path[0...-1]
              #log.debug "Enqueuing #{new_trail.to_shorthand} and [#{new_second_half.collect{|o| o.to_shorthand}.join(',') }]" if log.debug?
              problems_to_yield.push [new_trail, new_second_half]
            end
          end
        end
      end
    end

    def paths
      to_return = []
      each_path do |path|
        to_return.push path
      end
      to_return
    end

    def converged?
      !@converging_oriented_node_settable.nil?
    end

    # Return the OrientedNode that converges this bubble, behaviour
    # undefined if bubble is not converged
    def converging_oriented_node
      @known_problems[@converging_oriented_node_settable][0].path[-1]
    end

    def to_shorthand
      shorts = []
      if converged?
        shorts = paths.collect{|path| path.to_shorthand}
      else
        @queue.each do |problem|
          shorts.push problem.to_shorthand
        end
      end
      return "{#{shorts.join('|') }}"
    end

    def reverse!
      @is_reverse ||= false
      @is_reverse = !@is_reverse
    end

    # This doesn't make sense unless this is a converged bubble and the index == -1
    # because otherwise there is multiple answers
    def [](index)
      raise unless index == -1
      return Bio::Velvet::Graph::OrientedNodeTrail::OrientedNode.new(
        @converging_oriented_node_settable[0],
        @converging_oriented_node_settable[1]
        )
    end

    # Return one trail that exemplifies the paths through this bubble.
    # Current method of path selection is simply greedy, taking the highest coverage node
    # at each fork (or failing that the node with the lower node_id).
    def reference_trail
      raise unless converged?
      trail = []

      comparator = lambda do |problem1, problem2|
        node1 = nil
        node2 = nil
        if problem1.path.length == 1 and problem2.path.length > 1
          # Here the comparison cannot be made on node coverages
          # since one of the paths goes straight from the initial to the terminal
          # nodes. Choose instead based on if the second last node has higher or lower
          # coverage than the final node
          node1 = problem1.path[-1].node
          node2 = problem2.path[-2].node
        elsif problem2.path.length == 1 and problem1.path.length > 1
          node1 = problem1.path[-2].node
          node2 = problem2.path[-1].node
        else
          node1 = problem1.path[-2].node
          node2 = problem2.path[-2].node
        end

        if node1.coverage == node2.coverage
          -(node1.node_id <=> node2.node_id)
        else
          node1.coverage <=> node2.coverage
        end
      end

      current_problem = @known_problems[@converging_oriented_node_settable].max do |problem1, problem2|
        comparator.call problem1, problem2
      end
      reference_trail = []

      while !current_problem.nil?
        onode_to_add = current_problem.path[-1]
        raise CircuitousPathDetected if reference_trail.include?(onode_to_add)
        reference_trail.push onode_to_add

        if current_problem.path.length == 1
          current_problem = nil
        else
          current_problem = @known_problems[current_problem.path[-2].to_settable].max do |problem1, problem2|
            comparator.call problem1, problem2
          end
        end
      end

      node_trail = Bio::Velvet::Graph::OrientedNodeTrail.new
      node_trail.trail = reference_trail.reverse
      node_trail = node_trail.reverse if @is_reverse

      return node_trail
    end

    # Length in base pairs of the reference_trail
    def length_in_bp
      reference_trail.length_in_bp
    end

    # Does this (coverged) bubble contain any circuits?
    def circuitous?
      begin
        reference_trail #TODO: expand this to circuits within bubbles that aren't on the reference trail?
      rescue CircuitousPathDetected
        return true
      end
      return false
    end
  end

  class DynamicProgrammingProblem
    attr_accessor :path, :ubiquitous_oriented_nodes, :distance

    def initialize
      @path = []
      @ubiquitous_oriented_nodes = Set.new
    end

    def to_settable
      @path[-1].to_settable
    end

    def to_s
      ubiquitous_nodes = @ubiquitous_oriented_nodes.collect do |settabled|
        "#{settabled[0] }#{settabled[1] == Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST ? 's' : 'e'}"
      end
      return "DPP #{self.object_id}: #{@path.to_shorthand}/#{ubiquitous_nodes.join(',') }/#{distance}"
    end

    def to_shorthand
      ubiquitous_nodes = @ubiquitous_oriented_nodes.collect do |settabled|
        "#{settabled[0] }#{settabled[1] == Bio::Velvet::Graph::OrientedNodeTrail::START_IS_FIRST ? 's' : 'e'}"
      end
      "#{@path.to_shorthand}/#{ubiquitous_nodes.join(',') }/#{distance}"
    end
  end

  class CircuitousPathDetected < Exception; end
end
