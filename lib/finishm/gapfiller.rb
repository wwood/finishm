class Bio::FinishM::GapFiller
  def add_options(optparse_object, options)
    optparse_object.banner = "
    finishm gapfill <options..>
    "
  end

  def validate_options(options, argv)
    #return nil
    return 'I dunno'
  end
end
